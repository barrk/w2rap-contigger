///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2013) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "Set.h"
#include "VecUtilities.h"
#include "efasta/EfastaTools.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/HyperEfasta.h"
#include "paths/long/CleanEfasta.h"
#include "paths/long/Heuristics.h"
#include "paths/long/Logging.h"
#include "paths/long/LongProtoTools.h"
#include "paths/long/SupportedHyperBasevector.h"
#include "paths/long/Variants.h"




template<class H> void FlagEdgesForHiding( 
     const H& he, const vec<int>& inv, vec<Bool>& hide, const long_logging& logc )
{    double clock = WallClockTime( );
     vec<int> to_left, to_right;
     he.ToLeft(to_left), he.ToRight(to_right);

     // We break edges into five categories:
     // - keep_auto: not mapped by involution, or palindrome
     // - discard: to be hidden
     // - to_decide: edges not yet processed
     // - keep_look: edges to keep, but whose neighbors have not yet been looked at
     // - keep_looked: edges to keep whose neighbors have been examined.
     
     std::set<int> keep_auto, to_decide, keep_looked, keep_look, discard;

     // Define keep_auto, and initial to_decide.

     for ( int e = 0; e < he.EdgeObjectCount( ); e++ )
     {    if ( inv[e] < 0 || inv[e] == e ) keep_auto.insert(e);
          else to_decide.insert(e);    }

     // Main loop.

     while( keep_look.size( ) > 0 || to_decide.size( ) > 0 )
     {    if ( keep_look.size( ) > 0 )
          {    int e = *keep_look.begin( );
               keep_look.erase(e);
               keep_looked.insert(e);
               vec<int> nhood;
               int v = to_left[e], w = to_right[e];
               nhood.append( he.FromEdgeObj(v) );
               nhood.append( he.ToEdgeObj(v) );
               nhood.append( he.FromEdgeObj(w) );
               nhood.append( he.ToEdgeObj(w) );
               UniqueSort(nhood);
               for ( int j = 0; j < nhood.isize( ); j++ )
               {    int n = nhood[j];
                    if ( !Member( to_decide, n ) ) continue;
                    keep_look.insert(n), discard.insert( inv[n] );
                    to_decide.erase(n), to_decide.erase( inv[n] );    }    }
          else
          {    int e = *to_decide.begin( );
               to_decide.erase(e), to_decide.erase( inv[e] );
               keep_look.insert(e), discard.insert( inv[e] );    }    }

     // Report results.

     hide.resize( he.EdgeObjectCount( ), False );
     for ( std::set<int>::iterator i = discard.begin( ); i != discard.end( ); i++ )
          hide[*i] = True;
     REPORT_TIME( clock, "used flagging edges for hiding" );    }

template void FlagEdgesForHiding(
     const HyperEfasta& he, const vec<int>& inv, vec<Bool>& hide, 
     const long_logging& logc );

template void FlagEdgesForHiding(
     const SupportedHyperBasevector& he, const vec<int>& inv, vec<Bool>& hide,
     const long_logging& logc );

void CollapseBubbles( HyperEfasta& he )
{    const int max_brackets = 4;
     while(1)
     {    Bool changed = False;
          for ( int v = 0; v < he.N( ); v++ )
          {    for ( int l = 0; l < he.From(v).isize( ); l++ )
               {    int w = he.From(v)[l];
                    if ( w == v ) continue;
                    vec<int> ms, es;
                    for ( int m = 0; m < he.From(v).isize( ); m++ )
                    {    if ( he.From(v)[m] == w ) 
                         {    ms.push_back(m);
                              es.push_back( 
                                   he.EdgeObjectIndexByIndexFrom( v, m ) );    }    }
                    if ( es.solo( ) ) continue;

                    // Now we have two or more edges from vertex v to vertex w.

                    VecEFasta x;
                    for ( int m = 0; m < es.isize( ); m++ )
                         x.push_back( he.EdgeObject( es[m] ) );
                    int brackets = 0;
                    for ( size_t j = 0; j != x.size( ); j++ )
                    for ( int i = 0; i < x[j].isize( ); i++ )
                         if ( x[j][i] == '{' ) brackets++;
                    if ( brackets > max_brackets ) continue;
                    efasta e(x);
                    for ( int r = ms.isize( ) - 1; r >= 1; r-- )
                         he.DeleteEdgeFrom( v, ms[r] );
                    he.EdgeObjectByIndexFromMutable( v, ms[0] ) = e;
                    changed = True;
                    break;    }    }
          if ( !changed ) break;
          he.RemoveUnneededVertices( );    }
     he.RemoveDeadEdgeObjects( );
     he.RemoveEdgelessVertices( );    }

void GetCells( const HyperEfasta& he, vec<vec<int>>& cells )
{    cells.clear( );

     // Go through the vertices.

     for ( int v = 0; v < he.N( ); v++ )
     {    if ( he.From(v).solo( ) ) continue;

          // Now we have a vertex v, with at least two edges emanating from it.

          vec<int> suc, pre, between, betweenp;
          he.GetSuccessors1( v, suc );
          suc.EraseValue(v);
          for ( int wi = 0; wi < suc.isize( ); wi++ )
          {    int w = suc[wi];
               if ( he.To(w).solo( ) ) continue;

               // Now we have a vertex w which is not v, but follows v.  Moreover w
               // has at least two edges leading into it.  We posit v and w as
               // (possibly) bounding a cell.

               he.GetPredecessors1( w, pre );
               pre.EraseValue(w);
               between = Intersection( suc, pre );

               // Now between is the set of all vertices between v and w,
               // exclusive of them.

               betweenp = between;
               betweenp.push_back( v, w );
               Sort(betweenp);
               Bool bad = False;
               for ( int j = 0; j < between.isize( ); j++ )
               {    int x = between[j];
                    if ( !BinSubset( he.To(x), betweenp ) ) bad = True;
                    if ( !BinSubset( he.From(x), betweenp ) ) bad = True;    }
               if ( !BinSubset( he.From(v), betweenp ) ) bad = True;
               if ( !BinSubset( he.To(w), betweenp ) ) bad = True;
               if (bad) continue;

               // Now we know that v and w bound a cell.

               vec<int> cell;
               cell.push_back(v);
               cell.append(between);
               cell.push_back(w);
               cells.push_back(cell);    }    }    }

void Reduce( HyperEfasta& he, const int verbosity, const long_logging& logc )
{
     // First collapse bubbles.

     if ( verbosity >= 1 )
     {    std::cout << Date( ) << ": have " << he.EdgeObjectCount( ) << " edges, "
               << "collapsing bubbles" << std::endl;    }
     CollapseBubbles(he);

     // Now collapse more complex cells.

     if ( verbosity >= 1 ) 
     {    std::cout << Date( ) << ": have " << he.EdgeObjectCount( ) << " edges, "
               << "looking for cells to collapse" << std::endl;    }
     const int max_paths = 50;
     while(1)
     {    Bool progress = False;
          vec< vec<int> > cells;

          // Print.

          if ( verbosity >= 3 )
          {    std::cout << "\nat top of Reduce loop, graph is:\n";
               for ( int v = 0; v < he.N( ); v++ )
               for ( int j = 0; j < he.From(v).isize( ); j++ )
               {    int w = he.From(v)[j];
                    int e = he.EdgeObjectIndexByIndexFrom( v, j );
                    std::cout << v << " --(" << e << ",kmers=" 
                         << he.EdgeLengthKmers(e) 
                         << ", amb_events = " << he.EdgeObject(e).AmbEventCount( )
                         << ")--> " << w << std::endl;    }
               std::cout << std::endl;    }

          // Get cells.

          GetCells( he, cells );
          if ( verbosity >= 1 ) 
               std::cout << Date( ) << ": found " << cells.size( ) << " cells" << std::endl;

          // Go through the cells.

          for ( int i = 0; i < cells.isize( ); i++ )
          {    const vec<int>& cell = cells[i];
               if ( verbosity >= 2 )
                    std::cout << "examining cell " << printSeq(cell) << std::endl;
               if ( verbosity >= 3 )
               {    std::cout << "\ngraph is:\n";
                    for ( int v = 0; v < he.N( ); v++ )
                    for ( int j = 0; j < he.From(v).isize( ); j++ )
                    {    int w = he.From(v)[j];
                         int e = he.EdgeObjectIndexByIndexFrom( v, j );
                         std::cout << v << " --(" << e << ",kmers=" 
                              << he.EdgeLengthKmers(e) 
                              << ", amb_events = " 
                              << he.EdgeObject(e).AmbEventCount( )
                              << ")--> " << w << std::endl;    }
                    std::cout << std::endl;    }

               // Test for loop at w.  This is in effect working around a bug in
               // EdgePaths, which incorrectly ignores loops at w.

               if ( he.LoopAt( cell.back( ) ) ) continue;

               vec< vec<int> > paths;
               const int max_iterations = 1000;
               // note could upgrade to faster version that uses to_left, to_right
               Bool OK = he.EdgePaths( cell.front( ), cell.back( ), paths, 
                         -1, max_paths, max_iterations );
               if ( verbosity >= 1 ) PRINT2( int(OK), paths.size( ) );
               if ( !OK ) continue;
               if ( paths.isize( ) > max_paths ) continue;
               if ( paths.size( ) <= 1 ) continue;
               Bool bad = False;
               vec<int> edges;
               for ( int j = 0; j < paths.isize( ); j++ )
               for ( int k = 0; k < paths[j].isize( ); k++ )
                    edges.push_back( paths[j][k] );
               UniqueSort(edges);
               if ( verbosity >= 2 )
                    std::cout << "see edges " << printSeq(edges) << std::endl;
               VecEFasta bpaths;
               Bool fail = False;
               int K = he.K( );
               for ( int j = 0; j < paths.isize( ); j++ )
               {    if (fail) break;
                    efasta b = he.EdgeObject( paths[j][0] );
                    for ( int k = 1; k < paths[j].isize( ); k++ )
                    {    efasta c = he.EdgeObject( paths[j][k] );
                         if ( b.isize( ) < K - 1 )
                         {    fail = True;
                              break;    }
                         for ( int i = 0; i < K - 1; i++ )
                         {    if ( b[ b.isize( ) - 1 - i ] == '}' )
                              {    fail = True;
                                   break;    }    }
                         if ( !fail ) b.resize( b.isize( ) - ( K - 1 ) );
                         else
                         {    fail = False;
                              if ( c.isize( ) >= K )
                              {    for ( int i = 0; i < K - 1; i++ )
                                        if ( c[i] == '{' ) fail = True;    }
                              if ( !fail ) 
                                   c = c.substr( K-1, c.isize( ) - (K-1) );    }
                         if (fail) break;
                         b += c;    }
                    bpaths.push_back(b);    }
               if (fail) continue;
               int bracks = 0;
               for ( size_t i = 0; i != bpaths.size( ); i++ )
               for ( int j = 0; j < bpaths[i].isize( ); j++ )
                    if ( bpaths[i][j] == '{' ) bracks++;
               const int max_bracks = 10;
               if ( bracks > max_bracks ) continue;
               if ( verbosity >= 2 ) std::cout << "combining" << std::endl;
               efasta e(bpaths);
               he.DeleteEdges(edges);
               he.AddEdge( cell.front( ), cell.back( ), e );
               progress = True;
               he.RemoveUnneededVertices( );
               he.RemoveDeadEdgeObjects( );
               he.RemoveEdgelessVertices( );
               break;    }
          if ( !progress ) break;    }

     // Collapse bubbles again.

     CollapseBubbles(he);    

     // Find small unresolved but isolated cyclic edge clumps and 'box' them up in
     // a fastg-like contruction (no longer efasta).  Currently very restrictive.
     //
     // First define heuristics.

     const int min_long = 1000;
     const int max_short = 450;
     const int max_dist = 4;
     const int max_int = 15;

     // Repeat until no further progress.

     while(1)
     {    Bool improved = False;
          Bool deflower_verbose = ( logc.verb[ "DEFLOWER" ] >= 1 );
          if (deflower_verbose)
          {    std::cout << "\n";
               for ( int v = 0; v < he.N( ); v++ )
               for ( int j = 0; j < he.From(v).isize( ); j++ )
               {    int w = he.From(v)[j];
                    int e = he.EdgeObjectIndexByIndexFrom( v, j );
                    int len = he.EdgeLengthKmers(e);
                    PRINT4( v, w, e, len );    }
               std::cout << "\n";    }

          // Look for a long edge that will function as a left boundary.

          int K = he.K( );
          for ( int u = 0; u < he.N( ); u++ )
          for ( int j = 0; j < he.From(u).isize( ); j++ )
          {    int v = he.From(u)[j];
               int e1 = he.EdgeObjectIndexByIndexFrom( u, j );
               int n1 = he.EdgeObject(e1).isize( ) - K + 1;
               if ( n1 < min_long ) continue;
               if (deflower_verbose) PRINT(e1);
               
               // The left edge cannot have a bracket near its right end.
     
               Bool bad = False;
               for ( int s = 0; s < K - 1; s++ )
               {    if ( he.EdgeObject(e1)[ 
                         he.EdgeObject(e1).isize( ) - s - 1 ] == '}' )
                    {    bad = True;     }    }
               if (bad) continue;

               // Look for all vertices within distance d of v.

               vec< std::pair<int,int> > vd;
               vd.push( v, 0 );
               while(1)
               {    Bool progress = False;
                    for ( int i = 0; i < vd.isize( ); i++ )
                    {    int x = vd[i].first, d = vd[i].second;
                         if ( d == max_dist ) continue;
                         for ( int l = 0; l < he.From(x).isize( ); l++ )
                         {    int y = he.From(x)[l];
                              Bool found = False;
                              for ( int m = 0; m < vd.isize( ); m++ )
                              {    if ( vd[m].first == y )
                                   {    if ( d + 1 < vd[m].second )
                                        {    vd[m].second = d + 1;  
                                             progress = True;    }
                                        found = True;    }    }
                              if ( !found )
                              {    vd.push( y, d + 1 );
                                   progress = True;    }    }    }
                    if ( !progress ) break;    }
     
               // Now look for a long edge that could function as a right boundary.
     
               Bool edited = False;
               for ( int i = 0; i < vd.isize( ); i++ )
               {    if (edited) break;
                    int w = vd[i].first;
                    for ( int l = 0; l < he.From(w).isize( ); l++ )
                    {    int e2 = he.EdgeObjectIndexByIndexFrom( w, l );
                         if ( e2 == e1 ) continue;
                         int n2 = he.EdgeObject(e2).isize( ) - K + 1;
                         if ( n2 < min_long ) continue;
                         if (deflower_verbose) PRINT(e2);
     
                         // Now we have long edges that might function as left/right 
                         // boundaries.  Define what's between them.
          
                         vec<int> vint, eint;
                         vint.push_back( v, w );
                         //Bool get_back = False;
                         for ( int m = 0; m < vint.isize( ); m++ )
                         {    int x = vint[m];
                              for ( int r = 0; r < he.From(x).isize( ); r++ )
                              {    int y = he.From(x)[r];
                                   //if ( y == v ) get_back = True;
                                   int e = he.EdgeObjectIndexByIndexFrom( x, r );
     
                                   // Don't follow the right edge.
     
                                   if ( e == e2 ) continue;
     
                                   // Save vertex and edge.
     
                                   if (deflower_verbose)
                                        std::cout << "using edge " << e << std::endl;
                                   eint.push_back(e);
                                   if ( !Member( vint, y ) ) vint.push_back(y);
                                   if ( vint.isize( ) > max_int ) break;    }    }
                         if (deflower_verbose) PRINT( vint.size( ) );
     
                         // We have to be able to get back to v, and there can't be
                         // too many vertices.
     
                         if ( vint.isize( ) > max_int /* || !get_back */ ) continue;
     
                         // Going backwards can't take us out of vint.  Keep adding 
                         // to edge list.
     
                         Bool bad = False;
                         for ( int m = 0; m < vint.isize( ); m++ )
                         {    int x = vint[m];
                              for ( int r = 0; r < he.To(x).isize( ); r++ )
                              {    int e = he.EdgeObjectIndexByIndexTo( x, r );
                                   if ( e == e1 ) continue;
                                   eint.push_back(e);
                                   int y = he.To(x)[r];
                                   if ( !Member( vint, y ) ) bad = True;    }    }
                         if (deflower_verbose) PRINT(int(bad));
                         if (bad) continue;
     
                         // Make sure that edges are short.
     
                         UniqueSort(eint);
                         if (deflower_verbose) 
                              std::cout << "eint = {" << printSeq(eint) << "}" << std::endl;
                         Bool OK = True;
                         for ( int r = 0; r < eint.isize( ); r++ )
                         {    if ( he.EdgeLengthKmers( eint[r] ) > max_short ) 
                                   OK = False;    }
                         if (deflower_verbose) PRINT(int(OK));
                         if ( !OK ) continue;

                         // Create pseudo-fastg for the intermediate stuff.

                         if (deflower_verbose) std::cout << "deflowering" << std::endl;
                         if ( vint.size( ) > 2 )
                              std::swap( vint[1], vint[ vint.isize( ) - 1 ] );
                         String mid = "{";
                         for ( int m = 0; m < vint.isize( ); m++ )
                         {    int x = vint[m];
                              for ( int r = 0; r < he.From(x).isize( ); r++ )
                              {    int y = he.From(x)[r];
                                   int e = he.EdgeObjectIndexByIndexFrom( x, r );
                                   if ( e == e2 ) continue;
                                   int n = Position( vint, y );
                                   vec<basevector> b;
                                   he.EdgeObject(e).ExpandTo(b);
                                   for ( int z = 0; z < b.isize( ); z++ )
                                   {    if ( mid.size( ) > 1 ) mid += ",";
                                        mid += "(" + ToString(m) + "," 
                                        + ToString(n) + ",";
                                        for ( int s = 0; 
                                             s < b[z].isize( ) - ( K - 1 ); s++ )
                                        {    mid += as_base( b[z][s] );    }    
                                        mid += ")";
                                   }
                              }
                         }
                         mid += "}";

                         // Edit the assembly.
     
                         String left = he.EdgeObject(e1), right = he.EdgeObject(e2);
                         left.resize( left.isize( ) - ( K - 1 ) );
                         he.AddEdge( u, he.From(w)[l], left + mid + right );
                         vec<int> dels = eint;
                         dels.push_back( e1, e2 );
                         he.DeleteEdges(dels);
                         edited = True;
                         improved = True;
                         break;    }    }    }
          if ( !improved ) break;    }

     // Clean up.

     he.RemoveDeadEdgeObjects( );
     he.RemoveEdgelessVertices( );    }
