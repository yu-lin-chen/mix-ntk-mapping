#ifndef CHOICE_TO_LUTS_HPP
#define CHOICE_TO_LUTS_HPP
#pragma once

#include <mockturtle/mockturtle.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/utils/node_map.hpp>
#include <optional>
#include <unordered_map>
#include <unordered_set>
using namespace mockturtle;

namespace lsmap {

template <class NtkDest, class Ntk>
class collapse_mapped_choice_impl {
 public:
  using node_t   = typename Ntk::node;
  using signal_t = typename NtkDest::signal;

  collapse_mapped_choice_impl(Ntk const& ntk) : ntk(ntk) {}

  void run(NtkDest& dest) {
    mockturtle::node_map<signal_t, Ntk> node_to_signal(ntk);

    enum class driver_type { none, pos, neg, mixed };
    mockturtle::node_map<driver_type, Ntk> node_driver_type(ntk, driver_type::none);

    std::unordered_map<node_t, signal_t> opposites;

    // ===== helper: 统一处理 node/signal =====
    auto to_node = [&](auto x) -> node_t {
      using T = decltype(x);
      if constexpr (std::is_same_v<T, node_t>) {
        return x;
      } else {
        return ntk.get_node(x);
      }
    };

    auto is_compl = [&](auto x) -> bool {
      using T = decltype(x);
      if constexpr (std::is_same_v<T, node_t>) {
        return false;
      } else {
        return ntk.is_complemented(x);
      }
    };

    // ===== debug root/leaves =====
    ntk.foreach_node([&](auto const& n) {
      if (ntk.is_constant(n) || ntk.is_pi(n) || !ntk.is_cell_root(n)) return;

      std::cout << "root " << ntk.node_to_index(n) << " leaves:";
      ntk.foreach_cell_fanin(n, [&](auto fanin) {
        auto leaf = to_node(fanin);
        std::cout << " " << ntk.node_to_index(leaf);
      });
      std::cout << std::endl;
    });

    // ===== driver_type from PO =====
    ntk.foreach_po([&](auto const& f) {
      switch (node_driver_type[f]) {
        case driver_type::none:
          node_driver_type[f] = is_compl(f) ? driver_type::neg : driver_type::pos;
          break;
        case driver_type::pos:
          node_driver_type[f] = is_compl(f) ? driver_type::mixed : driver_type::pos;
          break;
        case driver_type::neg:
          node_driver_type[f] = is_compl(f) ? driver_type::neg : driver_type::mixed;
          break;
        case driver_type::mixed:
        default:
          break;
      }
    });

    // ===== driver_type from fanin =====
    ntk.foreach_node([&](auto const n) {
      if (ntk.is_constant(n) || ntk.is_pi(n) || !ntk.is_cell_root(n)) return;

      ntk.foreach_cell_fanin(n, [&](auto fanin) {
        auto leaf = to_node(fanin);
        if (node_driver_type[leaf] == driver_type::neg) {
          node_driver_type[leaf] = driver_type::mixed;
        }
      });
    });

    // ===== constants =====
    auto add_constant_to_map = [&](bool value) {
      const auto n = ntk.get_node(ntk.get_constant(value));
      switch (node_driver_type[n]) {
        default:
        case driver_type::none:
        case driver_type::pos:
          node_to_signal[n] = dest.get_constant(value);
          break;

        case driver_type::neg:
          node_to_signal[n] = dest.get_constant(!value);
          break;

        case driver_type::mixed:
          node_to_signal[n] = dest.get_constant(value);
          opposites[n] = dest.get_constant(!value);
          break;
      }
    };

    add_constant_to_map(false);
    if (ntk.get_node(ntk.get_constant(false)) !=
        ntk.get_node(ntk.get_constant(true))) {
      add_constant_to_map(true);
    }

    // ===== PIs =====
    ntk.foreach_pi([&](auto n) {
      signal_t dest_signal;
      switch (node_driver_type[n]) {
        default:
        case driver_type::none:
        case driver_type::pos:
          dest_signal = dest.create_pi();
          node_to_signal[n] = dest_signal;
          break;

        case driver_type::neg:
          dest_signal = dest.create_pi();
          node_to_signal[n] = dest.create_not(dest_signal);
          break;

        case driver_type::mixed:
          dest_signal = dest.create_pi();
          node_to_signal[n] = dest_signal;
          opposites[n] = dest.create_not(node_to_signal[n]);
          break;
      }

      if constexpr (mockturtle::has_has_name_v<Ntk> &&
                    mockturtle::has_get_name_v<Ntk> &&
                    mockturtle::has_set_name_v<NtkDest>) {
        if (ntk.has_name(ntk.make_signal(n)))
          dest.set_name(dest_signal, ntk.get_name(ntk.make_signal(n)));
      }
    });

    // ===== mapping 依赖递归构建 LUT =====
    std::unordered_set<node_t> built;

    std::function<void(node_t)> build = [&](node_t n) {
      if (built.count(n)) return;
      if (ntk.is_constant(n) || ntk.is_pi(n) || !ntk.is_cell_root(n)) return;

      // 先建 leaf
      ntk.foreach_cell_fanin(n, [&](auto fanin) {
        auto leaf = to_node(fanin);
        if (ntk.is_cell_root(leaf)) {
          build(leaf);
        }
      });

      // 再建当前 LUT
      std::vector<signal_t> children;
      ntk.foreach_cell_fanin(n, [&](auto fanin) {
        auto leaf = to_node(fanin);
        auto sig = node_to_signal[leaf];
        if (is_compl(fanin)) {
          sig = dest.create_not(sig);
        }
        children.push_back(sig);
      });

      switch (node_driver_type[n]) {
        default:
        case driver_type::none:
        case driver_type::pos:
          node_to_signal[n] = dest.create_node(children, ntk.cell_function(n));
          break;
        case driver_type::neg:
          node_to_signal[n] = dest.create_node(children, ~ntk.cell_function(n));
          break;
        case driver_type::mixed:
          node_to_signal[n] = dest.create_node(children, ntk.cell_function(n));
          opposites[n] = dest.create_node(children, ~ntk.cell_function(n));
          break;
      }

      built.insert(n);
    };

    // 从 PO 出发构建所有可达 root
    ntk.foreach_po([&](auto const& f) {
      auto root = to_node(f);
      if (ntk.is_cell_root(root)) {
        build(root);
      }
    });

    // ===== 输出 PO =====
    ntk.foreach_po([&](auto const& f, auto index) {
      (void)index;

      if (is_compl(f) && node_driver_type[f] == driver_type::mixed) {
        dest.create_po(opposites[to_node(f)]);
      } else {
        dest.create_po(node_to_signal[f]);
      }

      if constexpr (mockturtle::has_has_output_name_v<Ntk> &&
                    mockturtle::has_get_output_name_v<Ntk> &&
                    mockturtle::has_set_output_name_v<NtkDest>) {
        if (ntk.has_output_name(index)) {
          dest.set_output_name(index, ntk.get_output_name(index));
        }
      }
    });
  }

 private:
  Ntk const& ntk;
};

template <class NtkDest, class Ntk>
std::optional<NtkDest> choice_to_luts(Ntk const& ntk) {
  if (!ntk.has_mapping() && ntk.num_gates() > 0) {
    return std::nullopt;
  } else {
    lsmap::collapse_mapped_choice_impl<NtkDest, Ntk> p(ntk);
    NtkDest dest;
    p.run(dest);
    return dest;
  }
}

template <class NtkDest, class Ntk>
bool choice_to_luts(NtkDest& dest, Ntk const& ntk) {
  if (!ntk.has_mapping() && ntk.num_gates() > 0) {
    return false;
  } else {
    lsmap::collapse_mapped_choice_impl<NtkDest, Ntk> p(ntk);
    p.run(dest);
    return true;
  }
}

}  // namespace lsmap
#endif