#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <stack>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <kitty/dynamic_truth_table.hpp>
#include <kitty/operators.hpp>
#include <kitty/constructors.hpp>
#include <mockturtle/networks/aig.hpp>
#include <mockturtle/networks/events.hpp>
#include <mockturtle/traits.hpp>
#include <mockturtle/utils/algorithm.hpp>

namespace mockturtle
{

class mix_network
{
public:
  static constexpr auto min_fanin_size = 2u;
  static constexpr auto max_fanin_size = 3u;

  using base_type = mix_network;
  using node = uint32_t;

  enum class gate_type : uint8_t
  {
    none = 0,   // const / PI
    and2,
    xor2,
    or2,
    maj3,
    xor3,
    ite3
  };

  static constexpr node MIX_NULL{ 0x7FFFFFFFFFFFFFFFu };

  struct signal
  {
    signal() = default;
    signal( uint64_t index, uint64_t complement )
        : complement( complement ), index( index ) {}
    explicit signal( uint64_t data ) : data( data ) {}

    union
    {
      struct
      {
        uint64_t complement : 1;
        uint64_t index : 63;
      };
      uint64_t data{0};
    };

    signal operator!() const { return signal( data ^ 1 ); }
    signal operator+() const { return { index, 0 }; }
    signal operator-() const { return { index, 1 }; }
    signal operator^( bool c ) const { return signal( data ^ ( c ? 1 : 0 ) ); }
    bool operator==( signal const& other ) const { return data == other.data; }
    bool operator!=( signal const& other ) const { return data != other.data; }
    bool operator<( signal const& other ) const { return data < other.data; }
  };

  struct node_type
  {
    gate_type type{gate_type::none};
    uint8_t fanin{0};
    std::array<signal, 3> children{};
    uint32_t fanout_dead{0}; // MSB as dead flag (compatible with xmg)
    uint32_t value{0};
    uint32_t visited{0};
    uint8_t is_pi{0};

    bool operator==( node_type const& o ) const
    {
      if ( type != o.type || fanin != o.fanin ) return false;
      for ( uint8_t i = 0; i < fanin; ++i )
        if ( children[i] != o.children[i] ) return false;
      return true;
    }
  };

  struct node_hash
  {
    size_t operator()( node_type const& n ) const noexcept
    {
      uint64_t h = static_cast<uint64_t>( n.type ) * 1315423911u + n.fanin;
      for ( uint8_t i = 0; i < n.fanin; ++i )
        h = h * 1315423911u + n.children[i].data;
      return std::hash<uint64_t>{}( h );
    }
  };

  /* ---- storage for mockturtle views ---- */
  struct mix_storage
  {
    std::vector<node_type> nodes;
    std::vector<node> inputs;
    std::vector<signal> outputs;
    std::unordered_map<node_type, node, node_hash> hash;
    std::vector<signal> base_roots;
    std::vector<std::vector<uint32_t>> choice_table;
    uint32_t trav_id{0};
  };

  using storage = std::shared_ptr<mix_storage>;

public:
  mix_network()
      : _storage( std::make_shared<mix_storage>() ),
        _events( std::make_shared<network_events<base_type>>() ),
        _nodes( _storage->nodes ),
        _inputs( _storage->inputs ),
        _outputs( _storage->outputs ),
        _hash( _storage->hash ),
        base_roots_( _storage->base_roots ),
        choice_table_( _storage->choice_table )
  {
    // node 0 = constant false
    _nodes.emplace_back();
  }

  explicit mix_network( aig_network const& aig )
      : _storage( std::make_shared<mix_storage>() ),
        _events( std::make_shared<network_events<base_type>>() ),
        _nodes( _storage->nodes ),
        _inputs( _storage->inputs ),
        _outputs( _storage->outputs ),
        _hash( _storage->hash ),
        base_roots_( _storage->base_roots ),
        choice_table_( _storage->choice_table )
  {
    _nodes.emplace_back();
    assign_from_aig( aig );
  }
mix_network( std::shared_ptr<mix_storage> storage )
      : _storage( storage ),
        _events( std::make_shared<network_events<base_type>>() ),
        _nodes( _storage->nodes ),
        _inputs( _storage->inputs ),
        _outputs( _storage->outputs ),
        _hash( _storage->hash ),
        base_roots_( _storage->base_roots ),
        choice_table_( _storage->choice_table )
  {
  }

  mix_network clone() const
  {
    return { std::make_shared<mix_storage>( *_storage ) };
  }

  void assign_from_aig( aig_network const& aig )
  {
    ensure_base_roots( aig.size() );

    std::vector<signal> node_to_signal( aig.size(), get_constant( false ) );

    node_to_signal[0] = get_constant( false );
    base_roots_[0]     = get_constant( false );

    aig.foreach_pi( [&]( auto const& n ) {
      auto s   = create_pi();
      auto idx = aig.node_to_index( n );
      node_to_signal[idx] = s;
      base_roots_[idx]    = s;
    } );

    aig.foreach_gate( [&]( auto const& n ) {
      std::vector<signal> fins;
      fins.reserve( 2 );
      aig.foreach_fanin( n, [&]( auto const& f ) {
        auto s = node_to_signal[aig.node_to_index( aig.get_node( f ) )];
        if ( aig.is_complemented( f ) ) s = !s;
        fins.push_back( s );
      } );
      assert( fins.size() == 2u );
      auto s   = create_and( fins[0], fins[1] );
      auto idx = aig.node_to_index( n );
      node_to_signal[idx] = s;
      base_roots_[idx]    = s;
    } );

    aig.foreach_po( [&]( auto const& f ) {
      auto s = node_to_signal[aig.node_to_index( aig.get_node( f ) )];
      if ( aig.is_complemented( f ) ) s = !s;
      create_po( s );
    } );
  }

  /* 永久 choice 表接口 */
  // void add_choice( uint32_t root_idx, uint32_t cand_idx )
  // {
  //   ensure_choice_slot( root_idx );
  //   auto& v = choice_table_[root_idx];
  //   if ( std::find( v.begin(), v.end(), cand_idx ) == v.end() )
  //     v.push_back( cand_idx );
  // }
void add_choice( uint32_t root_idx, uint32_t cand_idx )
{
  ensure_choice_slot( root_idx );
  ensure_choice_slot( cand_idx );

  auto& v = choice_table_[root_idx];
  if ( std::find( v.begin(), v.end(), cand_idx ) == v.end() )
    v.push_back( cand_idx );

  auto& vr = choice_table_[cand_idx];
  if ( std::find( vr.begin(), vr.end(), root_idx ) == vr.end() )
    vr.push_back( root_idx );
}
  template<class Emit>
  void foreach_choice( uint32_t root_idx, Emit&& emit ) const
  {
    if ( root_idx >= choice_table_.size() ) return;
    for ( auto cand : choice_table_[root_idx] ) emit( cand );
  }

  auto const& choice_table() const noexcept { return choice_table_; }

  bool is_repr( node const& n ) const
  {
    if ( n >= choice_table_.size() ) return false;
    return !choice_table_[n].empty() && fanout_size( n ) > 0;
  }

  node get_equiv_node( node const& n ) const
  {
    if ( n >= choice_table_.size() || choice_table_[n].empty() ) return MIX_NULL;
    return choice_table_[n].front();
  }
  // po--pool
void keep_outputs( std::vector<signal> const& keep )
{
  std::cout << "[mix] keep_outputs called: old=" << _outputs.size()
            << " new=" << keep.size() << "\n";

  std::cout << "  old outputs: ";
  for ( auto const& o : _outputs )
    std::cout << (o.complement ? "!" : "") << o.index << " ";
  std::cout << "\n";

  std::cout << "  new outputs: ";
  for ( auto const& o : keep )
    std::cout << (o.complement ? "!" : "") << o.index << " ";
  std::cout << "\n";

  for ( auto const& o : _outputs )
    decr_fanout_size( static_cast<node>( o.index ) );

  _outputs = keep;

  for ( auto const& o : _outputs )
    incr_fanout_size( static_cast<node>( o.index ) );
}
  // ===== constants / PI / PO =====
  signal get_constant( bool value ) const { return signal( 0, value ? 1 : 0 ); }

  signal create_pi()
  {
    node_type n;
    n.is_pi = 1;
    const auto idx = static_cast<node>( _nodes.size() );
    _nodes.push_back( n );
    _inputs.push_back( idx );
    return signal( idx, 0 );
  }

  uint32_t create_po( signal const& f )
  {
    _nodes[f.index].fanout_dead++;
    _outputs.push_back( f );
    return static_cast<uint32_t>( _outputs.size() - 1 );
  }

  bool is_combinational() const { return true; }

  bool constant_value( node const& ) const { return false; }

  // ===== unary =====
  signal create_buf( signal const& a ) { return a; }
  signal create_not( signal const& a ) { return !a; }

  // ===== create gates =====
  signal create_and( signal a, signal b )
  {
    if ( a.index > b.index ) std::swap( a, b );
    if ( a.index == b.index ) return ( a.complement == b.complement ) ? a : get_constant( false );
    if ( a.index == 0 ) return a.complement ? b : get_constant( false );
    return _create_node2( gate_type::and2, a, b );
  }

  signal create_or( signal a, signal b )
  {
    if ( a.index > b.index ) std::swap( a, b );
    if ( a.index == b.index ) return ( a.complement == b.complement ) ? a : get_constant( true );
    if ( a.index == 0 ) return a.complement ? get_constant( true ) : b;
    return _create_node2( gate_type::or2, a, b );
  }

  signal create_xor( signal a, signal b )
  {
    if ( a.index < b.index ) std::swap( a, b );
    bool out_c = a.complement != b.complement;
    a.complement = b.complement = 0;
    if ( a.index == b.index ) return get_constant( out_c );
    if ( b.index == 0 ) return a ^ out_c;
    return _create_node2( gate_type::xor2, a, b ) ^ out_c;
  }

  signal create_maj( signal a, signal b, signal c )
  {
    _sort3_asc( a, b, c );
    if ( a.index == b.index ) return ( a.complement == b.complement ) ? a : c;
    if ( b.index == c.index ) return ( b.complement == c.complement ) ? b : a;
    return _create_node3( gate_type::maj3, a, b, c );
  }

  signal create_xor3( signal a, signal b, signal c )
  {
    _sort3_desc( a, b, c );
    bool out_c = ( a.complement != b.complement ) != c.complement;
    a.complement = b.complement = c.complement = 0;
    if ( a.index == b.index ) return c ^ out_c;
    if ( b.index == c.index ) return a ^ out_c;
    return _create_node3( gate_type::xor3, a, b, c ) ^ out_c;
  }

  signal create_ite( signal cond, signal t, signal e )
  {
    return _create_node3( gate_type::ite3, cond, t, e );
  }

  signal create_nand( signal const& a, signal const& b ) { return !create_and( a, b ); }
  signal create_nor( signal const& a, signal const& b ) { return !create_or( a, b ); }
  signal create_lt( signal const& a, signal const& b ) { return create_and( !a, b ); }
  signal create_le( signal const& a, signal const& b ) { return !create_and( a, !b ); }
  signal create_xnor( signal const& a, signal const& b ) { return create_xor( a, b ) ^ true; }

  // n-ary
  signal create_nary_and( std::vector<signal> const& fs )
  {
    return tree_reduce( fs.begin(), fs.end(), get_constant( true ), [this]( auto const& a, auto const& b ) { return create_and( a, b ); } );
  }

  signal create_nary_or( std::vector<signal> const& fs )
  {
    return tree_reduce( fs.begin(), fs.end(), get_constant( false ), [this]( auto const& a, auto const& b ) { return create_or( a, b ); } );
  }

  signal create_nary_xor( std::vector<signal> const& fs )
  {
    return ternary_tree_reduce( fs.begin(), fs.end(), get_constant( false ), [this]( auto const& a, auto const& b, auto const& c ) { return create_xor3( a, b, c ); } );
  }

  signal clone_node( mix_network const& other, node const& source, std::vector<signal> const& children )
  {
    if ( other.is_maj( source ) ) return create_maj( children[0], children[1], children[2] );
    if ( other.is_xor3( source ) ) return create_xor3( children[0], children[1], children[2] );
    if ( other.is_and( source ) ) return create_and( children[0], children[1] );
    if ( other.is_or( source ) )  return create_or( children[0], children[1] );
    if ( other.is_xor( source ) ) return create_xor( children[0], children[1] );
    if ( other.is_ite( source ) ) return create_ite( children[0], children[1], children[2] );
    return children[0];
  }

  // ===== basic query =====
  uint32_t size() const { return static_cast<uint32_t>( _nodes.size() ); }
  uint32_t num_pis() const { return static_cast<uint32_t>( _inputs.size() ); }
  uint32_t num_pos() const { return static_cast<uint32_t>( _outputs.size() ); }
  uint32_t num_gates() const { return static_cast<uint32_t>( _hash.size() ); }
  signal base_root(uint32_t idx) const { return base_roots_.at(idx); }
  bool is_constant( node n ) const { return n == 0; }
  bool is_pi( node n ) const { return n < _nodes.size() && _nodes[n].is_pi && n != 0; }
  bool is_ci( node n ) const { return is_pi( n ); }

  bool is_and( node n ) const { return _nodes[n].type == gate_type::and2; }
  bool is_or ( node n ) const { return _nodes[n].type == gate_type::or2;  }
  bool is_xor( node n ) const { return _nodes[n].type == gate_type::xor2; }
  bool is_maj( node n ) const { return _nodes[n].type == gate_type::maj3; }
  bool is_xor3(node n ) const { return _nodes[n].type == gate_type::xor3; }
  bool is_ite( node n ) const { return _nodes[n].type == gate_type::ite3; }

  uint32_t node_to_index( node const& n ) const { return static_cast<uint32_t>( n ); }
  node index_to_node( uint32_t index ) const { return static_cast<node>( index ); }
  node get_node( signal const& s ) const { return static_cast<node>( s.index ); }
  signal make_signal( node n ) const { return signal( n, 0 ); }
  bool is_complemented( signal const& s ) const { return s.complement; }

  // ===== fanout / value / visited =====
  uint32_t fanout_size( node n ) const { return _nodes[n].fanout_dead & 0x7fffffff; }
  uint32_t incr_fanout_size( node n ) { return _nodes[n].fanout_dead++ & 0x7fffffff; }
  uint32_t decr_fanout_size( node n ) { return --_nodes[n].fanout_dead & 0x7fffffff; }

  bool is_dead( node n ) const { return ( _nodes[n].fanout_dead >> 31 ) & 1; }

  void clear_values() { for ( auto& n : _nodes ) n.value = 0; }
  uint32_t value( node n ) const { return _nodes[n].value; }
  void set_value( node n, uint32_t v ) { _nodes[n].value = v; }

  void clear_visited() { for ( auto& n : _nodes ) n.visited = 0; }
  uint32_t visited( node n ) const { return _nodes[n].visited; }
  void set_visited( node n, uint32_t v ) { _nodes[n].visited = v; }

  uint32_t trav_id() const { return _storage->trav_id; }
  void incr_trav_id() { ++_storage->trav_id; }

  // ===== traversal =====
  template<typename Fn>
  void foreach_node( Fn&& fn ) const
  {
    for ( node n = 0; n < _nodes.size(); ++n )
      if ( !is_dead( n ) ) fn( n );
  }

  template<typename Fn>
void foreach_ci( Fn&& fn ) const
{
  for ( uint32_t i = 0; i < _inputs.size(); ++i )
  {
    if constexpr ( detail::is_callable_with_index_v<Fn, node, void> ||
                   detail::is_callable_with_index_v<Fn, node, bool> )
    {
      if constexpr ( detail::is_callable_with_index_v<Fn, node, bool> )
      {
        if ( !fn( _inputs[i], i ) ) return;
      }
      else
      {
        fn( _inputs[i], i );
      }
    }
    else
    {
      if constexpr ( detail::is_callable_without_index_v<Fn, node, bool> )
      {
        if ( !fn( _inputs[i] ) ) return;
      }
      else
      {
        fn( _inputs[i] );
      }
    }
  }
}

template<typename Fn>
void foreach_co( Fn&& fn ) const
{
  for ( uint32_t i = 0; i < _outputs.size(); ++i )
  {
    if constexpr ( detail::is_callable_with_index_v<Fn, signal, void> ||
                   detail::is_callable_with_index_v<Fn, signal, bool> )
    {
      if constexpr ( detail::is_callable_with_index_v<Fn, signal, bool> )
      {
        if ( !fn( _outputs[i], i ) ) return;
      }
      else
      {
        fn( _outputs[i], i );
      }
    }
    else
    {
      if constexpr ( detail::is_callable_without_index_v<Fn, signal, bool> )
      {
        if ( !fn( _outputs[i] ) ) return;
      }
      else
      {
        fn( _outputs[i] );
      }
    }
  }
}

template<typename Fn>
void foreach_pi( Fn&& fn ) const
{
  foreach_ci( std::forward<Fn>( fn ) );
}

template<typename Fn>
void foreach_po( Fn&& fn ) const
{
  foreach_co( std::forward<Fn>( fn ) );
}

  template<typename Fn>
  void foreach_gate( Fn&& fn ) const
  {
    for ( node n = 1; n < _nodes.size(); ++n )
      if ( !is_pi( n ) && !is_dead( n ) ) fn( n );
  }

template<typename Fn>
void foreach_fanin( node n, Fn&& fn ) const
{
  if ( is_constant( n ) || is_pi( n ) ) return;

  static_assert( detail::is_callable_without_index_v<Fn, signal, bool> ||
                 detail::is_callable_with_index_v<Fn, signal, bool> ||
                 detail::is_callable_without_index_v<Fn, signal, void> ||
                 detail::is_callable_with_index_v<Fn, signal, void> );

  const auto k = fanin_size( n );
  for ( uint32_t i = 0; i < k; ++i )
  {
    if constexpr ( detail::is_callable_with_index_v<Fn, signal, bool> )
    {
      if ( !fn( _nodes[n].children[i], i ) ) return;
    }
    else if constexpr ( detail::is_callable_with_index_v<Fn, signal, void> )
    {
      fn( _nodes[n].children[i], i );
    }
    else if constexpr ( detail::is_callable_without_index_v<Fn, signal, bool> )
    {
      if ( !fn( _nodes[n].children[i] ) ) return;
    }
    else
    {
      fn( _nodes[n].children[i] );
    }
  }
}

  signal get_child0( node const& p ) const { return _nodes[p].children[0]; }
  signal get_child1( node const& p ) const { return _nodes[p].children[1]; }
  signal get_child2( node const& p ) const { return _nodes[p].children[2]; }

  uint32_t fanin_size( node n ) const
  {
    switch ( _nodes[n].type )
    {
    case gate_type::and2:
    case gate_type::or2:
    case gate_type::xor2: return 2;
    case gate_type::maj3:
    case gate_type::xor3:
    case gate_type::ite3: return 3;
    default: return 0;
    }
  }

  void replace_in_outputs( node old_n, signal const& new_s )
  {
    for ( auto& o : _outputs )
      if ( o.index == old_n )
      {
        o.index = new_s.index;
        o.complement ^= new_s.complement;
        if ( old_n != new_s.index ) incr_fanout_size( static_cast<node>( new_s.index ) );
      }
  }
  gate_type get_gate_type( node n ) const { return _nodes[n].type; }
  // ===== index helpers (align with xmg) =====
  node ci_at( uint32_t index ) const { return _inputs.at( index ); }
  signal co_at( uint32_t index ) const { return _outputs.at( index ); }
  node pi_at( uint32_t index ) const { return _inputs.at( index ); }
  signal po_at( uint32_t index ) const { return _outputs.at( index ); }

  uint32_t ci_index( node const& n ) const
  {
    for ( uint32_t i = 0; i < _inputs.size(); ++i ) if ( _inputs[i] == n ) return i;
    return static_cast<uint32_t>( -1 );
  }

  uint32_t co_index( signal const& s ) const
  {
    for ( uint32_t i = 0; i < _outputs.size(); ++i ) if ( _outputs[i] == s ) return i;
    return static_cast<uint32_t>( -1 );
  }

  uint32_t pi_index( node const& n ) const { return ci_index( n ); }
  uint32_t po_index( signal const& s ) const { return co_index( s ); }

  // ===== simulation (align xmg) =====
kitty::dynamic_truth_table node_function( const node& n ) const
{
  switch ( _nodes[n].type )
  {
  case gate_type::and2:
  {
    kitty::dynamic_truth_table tt( 2 );
    kitty::create_nth_var( tt, 0 ); // a
    auto b = tt; kitty::create_nth_var( b, 1 );
    return tt & b;
  }
  case gate_type::or2:
  {
    kitty::dynamic_truth_table tt( 2 );
    kitty::create_nth_var( tt, 0 );
    auto b = tt; kitty::create_nth_var( b, 1 );
    return tt | b;
  }
  case gate_type::xor2:
  {
    kitty::dynamic_truth_table tt( 2 );
    kitty::create_nth_var( tt, 0 );
    auto b = tt; kitty::create_nth_var( b, 1 );
    return tt ^ b;
  }
  case gate_type::maj3:
  {
    kitty::dynamic_truth_table a( 3 ), b( 3 ), c( 3 );
    kitty::create_nth_var( a, 0 );
    kitty::create_nth_var( b, 1 );
    kitty::create_nth_var( c, 2 );
    return kitty::ternary_majority( a, b, c );
  }
  case gate_type::xor3:
  {
    kitty::dynamic_truth_table a( 3 ), b( 3 ), c( 3 );
    kitty::create_nth_var( a, 0 );
    kitty::create_nth_var( b, 1 );
    kitty::create_nth_var( c, 2 );
    return a ^ b ^ c;
  }
  case gate_type::ite3:
  {
    kitty::dynamic_truth_table s( 3 ), t( 3 ), e( 3 );
    kitty::create_nth_var( s, 0 );
    kitty::create_nth_var( t, 1 );
    kitty::create_nth_var( e, 2 );
    return ( s & t ) | ( ~s & e );
  }
  default:
  {
    kitty::dynamic_truth_table tt( 1 );
    tt._bits[0] = 0;
    return tt;
  }
  }
}

template<typename Iterator>
iterates_over_t<Iterator, bool>
compute( node const& n, Iterator begin, Iterator end ) const
{
  (void)end;
  assert( n != 0 && !is_ci( n ) );

  auto const& c1 = _nodes[n].children[0];
  auto const& c2 = _nodes[n].children[1];

  auto v1 = *begin++;
  auto v2 = *begin++;

  switch ( _nodes[n].type )
  {
  case gate_type::and2:
    return ( v1 ^ c1.complement ) & ( v2 ^ c2.complement );
  case gate_type::or2:
    return ( v1 ^ c1.complement ) | ( v2 ^ c2.complement );
  case gate_type::xor2:
    return ( v1 ^ c1.complement ) ^ ( v2 ^ c2.complement );

  case gate_type::maj3:
  case gate_type::xor3:
  case gate_type::ite3:
  {
    auto const& c3 = _nodes[n].children[2];
    auto v3 = *begin++;
    if ( is_xor3( n ) )
      return ( ( v1 ^ c1.complement ) != ( v2 ^ c2.complement ) ) !=
             ( v3 ^ c3.complement );
    // maj3 / ite3 均用 maj 语义（保持原逻辑）
    return ( ( v1 ^ c1.complement ) && ( v2 ^ c2.complement ) ) ||
           ( ( v3 ^ c3.complement ) && ( v1 ^ c1.complement ) ) ||
           ( ( v3 ^ c3.complement ) && ( v2 ^ c2.complement ) );
  }
  default:
    break;
  }
  return false;
}

template<typename Iterator>
iterates_over_truth_table_t<Iterator>
compute( node const& n, Iterator begin, Iterator end ) const
{
  (void)end;
  assert( n != 0 && !is_ci( n ) );

  auto const& c1 = _nodes[n].children[0];
  auto const& c2 = _nodes[n].children[1];

  auto tt1 = *begin++;
  auto tt2 = *begin++;

  switch ( _nodes[n].type )
  {
  case gate_type::and2:
    return ( c1.complement ? ~tt1 : tt1 ) & ( c2.complement ? ~tt2 : tt2 );
  case gate_type::or2:
    return ( c1.complement ? ~tt1 : tt1 ) | ( c2.complement ? ~tt2 : tt2 );
  case gate_type::xor2:
    return ( c1.complement ? ~tt1 : tt1 ) ^ ( c2.complement ? ~tt2 : tt2 );

  case gate_type::maj3:
  case gate_type::xor3:
  case gate_type::ite3:
  {
    auto const& c3 = _nodes[n].children[2];
    auto tt3 = *begin++;
    if ( is_xor3( n ) )
      return ( c1.complement ? ~tt1 : tt1 ) ^
             ( c2.complement ? ~tt2 : tt2 ) ^
             ( c3.complement ? ~tt3 : tt3 );
    return kitty::ternary_majority( c1.complement ? ~tt1 : tt1,
                                    c2.complement ? ~tt2 : tt2,
                                    c3.complement ? ~tt3 : tt3 );
  }
  default:
    break;
  }
  return tt1; // fallback，不应到达
}

template<typename Iterator>
void compute( node const& n, kitty::partial_truth_table& result, Iterator begin, Iterator end ) const
{
  static_assert( iterates_over_v<Iterator, kitty::partial_truth_table>,
                 "begin and end have to iterate over partial_truth_tables" );
  (void)end;

  assert( n != 0 && !is_ci( n ) );

  auto const& c1 = _nodes[n].children[0];
  auto const& c2 = _nodes[n].children[1];

  auto tt1 = *begin++;
  auto tt2 = *begin++;

  result.resize( tt1.num_bits() );

  switch ( _nodes[n].type )
  {
  case gate_type::and2:
    result._bits.back() =
        ( c1.complement ? ~tt1._bits.back() : tt1._bits.back() ) &
        ( c2.complement ? ~tt2._bits.back() : tt2._bits.back() );
    break;

  case gate_type::or2:
    result._bits.back() =
        ( c1.complement ? ~tt1._bits.back() : tt1._bits.back() ) |
        ( c2.complement ? ~tt2._bits.back() : tt2._bits.back() );
    break;

  case gate_type::xor2:
    result._bits.back() =
        ( c1.complement ? ~tt1._bits.back() : tt1._bits.back() ) ^
        ( c2.complement ? ~tt2._bits.back() : tt2._bits.back() );
    break;

  case gate_type::maj3:
  case gate_type::xor3:
  case gate_type::ite3:
  {
    auto const& c3 = _nodes[n].children[2];
    auto tt3 = *begin++;
    if ( is_xor3( n ) )
    {
      result._bits.back() =
          ( c1.complement ? ~tt1._bits.back() : tt1._bits.back() ) ^
          ( c2.complement ? ~tt2._bits.back() : tt2._bits.back() ) ^
          ( c3.complement ? ~tt3._bits.back() : tt3._bits.back() );
    }
    else
    {
      result._bits.back() =
          ( ( c1.complement ? ~tt1._bits.back() : tt1._bits.back() ) &
            ( c2.complement ? ~tt2._bits.back() : tt2._bits.back() ) ) |
          ( ( c1.complement ? ~tt1._bits.back() : tt1._bits.back() ) &
            ( c3.complement ? ~tt3._bits.back() : tt3._bits.back() ) ) |
          ( ( c2.complement ? ~tt2._bits.back() : tt2._bits.back() ) &
            ( c3.complement ? ~tt3._bits.back() : tt3._bits.back() ) );
    }
    break;
  }
  default:
    break;
  }

  result.mask_bits();
}

  auto& events() const { return *_events; }

private:
  signal _create_node2( gate_type t, signal a, signal b )
  {
    node_type n;
    n.type = t;
    n.fanin = 2;
    n.children[0] = a;
    n.children[1] = b;

    if ( auto it = _hash.find( n ); it != _hash.end() )
      return signal( it->second, 0 );

    const auto idx = static_cast<node>( _nodes.size() );
    _nodes.push_back( n );
    _hash[_nodes.back()] = idx;

    incr_fanout_size( static_cast<node>( a.index ) );
    incr_fanout_size( static_cast<node>( b.index ) );
    return signal( idx, 0 );
  }

  signal _create_node3( gate_type t, signal a, signal b, signal c )
  {
    node_type n;
    n.type = t;
    n.fanin = 3;
    n.children[0] = a;
    n.children[1] = b;
    n.children[2] = c;

    if ( auto it = _hash.find( n ); it != _hash.end() )
      return signal( it->second, 0 );

    const auto idx = static_cast<node>( _nodes.size() );
    _nodes.push_back( n );
    _hash[_nodes.back()] = idx;

    incr_fanout_size( static_cast<node>( a.index ) );
    incr_fanout_size( static_cast<node>( b.index ) );
    incr_fanout_size( static_cast<node>( c.index ) );
    return signal( idx, 0 );
  }

  static void _sort3_asc( signal& a, signal& b, signal& c )
  {
    if ( a.index > b.index ) std::swap( a, b );
    if ( b.index > c.index ) std::swap( b, c );
    if ( a.index > b.index ) std::swap( a, b );
  }

  static void _sort3_desc( signal& a, signal& b, signal& c )
  {
    if ( a.index < b.index ) std::swap( a, b );
    if ( b.index < c.index ) std::swap( b, c );
    if ( a.index < b.index ) std::swap( a, b );
  }

private:
  std::shared_ptr<mix_storage> _storage;
  std::shared_ptr<network_events<base_type>> _events;

  // references to storage fields (keep old names / structure)
  std::vector<node_type>& _nodes;
  std::vector<node>& _inputs;
  std::vector<signal>& _outputs;
  std::unordered_map<node_type, node, node_hash>& _hash;
  std::vector<signal>& base_roots_;
  std::vector<std::vector<uint32_t>>& choice_table_;

  void ensure_base_roots(uint32_t sz)
  {
    if (base_roots_.size() < sz) base_roots_.resize(sz, get_constant(false));
  }

  void ensure_choice_slot(uint32_t idx)
  {
      if ( idx == MIX_NULL ) return;  // 防御
  if ( choice_table_.size() <= idx ) choice_table_.resize(idx + 1);
    if (choice_table_.size() <= idx) choice_table_.resize(idx + 1);
  }
};

} // namespace mockturtle

// ---- register network type ----
namespace mockturtle
{
template<>
struct is_network_type<mix_network> : std::true_type {};
}

namespace std
{
template<>
struct hash<mockturtle::mix_network::signal>
{
  uint64_t operator()( mockturtle::mix_network::signal const& s ) const noexcept
  {
    uint64_t k = s.data;
    k ^= k >> 33;
    k *= 0xff51afd7ed558ccd;
    k ^= k >> 33;
    k *= 0xc4ceb9fe1a85ec53;
    k ^= k >> 33;
    return k;
  }
};
} // namespace std