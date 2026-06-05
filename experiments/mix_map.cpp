#include <iostream>
#include <string>
#include <vector>
#include <fmt/format.h>
#include <lorina/aiger.hpp>
#include <mockturtle/mockturtle.hpp>
#include <../mixmap/choice_lut_mapper.hpp>
#include <../mixmap/choice_to_luts.hpp>
#include <fstream>
#include <mockturtle/io/write_blif.hpp>
using namespace mockturtle;

template<typename Ntk>
void dump_network_edges(const Ntk& ntk, const std::string& tag)
{
  std::cout << "\n=== " << tag << " ===\n";
  std::cout << fmt::format("pis={} pos={} gates={} size={}\n",
                           ntk.num_pis(), ntk.num_pos(), ntk.num_gates(), ntk.size());

  // constants
  auto c0 = ntk.get_node(ntk.get_constant(false));
  std::cout << fmt::format("const0: node {}\n", ntk.node_to_index(c0));

  // PIs
  std::cout << "PIs:\n";
  ntk.foreach_pi([&](auto n, auto i) {
    std::cout << fmt::format("  pi[{}] -> node {}\n", i, ntk.node_to_index(n));
  });

  // internal nodes (gates)
  std::cout << "Gates:\n";
  ntk.foreach_gate([&](auto n) {
    auto idx = ntk.node_to_index(n);
    std::vector<std::string> fins;
    ntk.foreach_fanin(n, [&](auto f, auto pin) {
      auto cn = ntk.get_node(f);
      auto cidx = ntk.node_to_index(cn);
      bool comp = ntk.is_complemented(f);
      fins.push_back(fmt::format("{}{}",
                                 comp ? "!" : "",
                                 cidx));
    });
    // std::cout << fmt::format("  node {} <- {}type {}\n", idx, fmt::join(fins, ", "),ntk.get_gate_type(n));
   std::cout << fmt::format("  node {} <- {}\n", idx, fmt::join(fins, ", "));
  });

  // POs
  std::cout << "POs:\n";
  ntk.foreach_po([&](auto f, auto i) {
    auto n = ntk.get_node(f);
    auto idx = ntk.node_to_index(n);
    bool comp = ntk.is_complemented(f);
    std::cout << fmt::format("  po[{}] <- {}{}\n", i, comp ? "!" : "", idx);
  });
}

int main(int argc, char** argv)
{
  if (argc < 2) {
    std::cerr << "usage: ./run <input.aig>\n";
    return 1;
  }
  std::string filename = argv[1];
  aig_network aig;
  if (lorina::read_aiger(filename, aiger_reader(aig)) != lorina::return_code::success) {
    std::cerr << "read_aiger failed: " << filename << "\n";
    return 1;
  }
  //aig-->mix 
  // dump_network_edges(aig, "AIG");
  mix_network mix{aig}; 
  //exact map AIG -> XMG  | create xmg --> mix
  exact_library_params eps;
  xmg_npn_resynthesis resyn;
  exact_library<xmg_network> lib(resyn, eps);
  map_params ps;
  map_stats st;
  xmg_network xmg = mockturtle::map(aig, lib, ps, &st,&mix);
    //exact map AIG -> MIG  | create mig --> mix
  // mig_npn_resynthesis resyn1;
  // exact_library<mig_network> lib1(resyn1, eps);
  // map_params ps1;
  // map_stats st1;
  // mig_network mig = mockturtle::map(aig, lib1, ps1, &st1,&mix);
    //exact map AIG -> XAG  | create xmg --> mix
  gate_dot_drawer<mix_network> drawer;
  // 写到文件
  std::ofstream os( "mix.dot" );
  write_dot( mix, os, drawer );
  os.close();
   mapping_view<mix_network, true> mapped{ mix };
   also::choice_lut_map_params pst;
pst.cross_choice_provider = [&mix](uint32_t idx, auto const& emit){
  // idx 是当前映射网络的节点 index，需要与 mix 的 base_roots_/节点 index 对齐
  if (idx >= mix.choice_table().size()) return;
  for (auto cand : mix.choice_table()[idx]) emit(cand);
};
 also::choice_lut_map_inplace(mapped, pst); 
 //mapped type:mapping_view mix_ntk
  klut_network klut = *lsmap::choice_to_luts<klut_network>( mapped );
//   klut.foreach_gate([&](auto n) {
//   std::cout << "LUT node " << n << ": ";
//   klut.foreach_fanin(n, [&](auto f) {
//     std::cout << klut.get_node(f) << " ";
//   });
//   std::cout << " func=" << kitty::to_hex(klut.node_function(n));
//   std::cout << std::endl;
// });
          mockturtle::depth_view depth_klut{ klut };
          std::cout << fmt::format( "LUT mapped XMG into #gates = {} level = {}\n",
                                    klut.num_gates(), depth_klut.depth() );
//   std::ofstream os("out.blif");
// write_blif( klut, os );

   // dump_network_edges(xmg, "XMG(mapped)");
  //     dump_network_edges(mig, "MIG(mapped)");
      dump_network_edges(mix, "mix");
    const auto& tbl = mix.choice_table();
  std::cout << "\n=== mix choice_table ===\n";
  for (uint32_t idx = 0; idx < tbl.size(); ++idx) {
    if (tbl[idx].empty()) continue;
    std::cout << "root " << idx << " :";
    for (auto cand : tbl[idx]) std::cout << " " << cand;
    std::cout << "\n";
  }
   aig_network aiger_k = convert_klut_to_graph<aig_network>(klut);
        const auto miter_xmg = *miter<xmg_network>( aig, aiger_k ); 
        equivalence_checking_stats eq_st;
        const auto result = equivalence_checking( miter_xmg, {}, &eq_st );
        assert( result );
        assert( *result );
        std::cout << "Network is equivalent after resub\n";
  return 0;
}
