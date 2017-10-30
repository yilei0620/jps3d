/**
 * @file graph_search.h
 * @brief backend of graph search, implemetation of A* and Lifelong Planning A*
 */

#ifndef GRAPH_SEARCH_H
#define GRAPH_SEARCH_H

#include <boost/heap/d_ary_heap.hpp>      // boost::heap::d_ary_heap
#include <memory>                         // std::shared_ptr
#include <limits>                         // std::numeric_limits
#include <vector>                         // std::vector
#include <unordered_map>                  // std::unordered_map
#include <array>                          // std::array

namespace JPS
{
  ///Declare `env' class
  class env_base;

  
  ///Heap element comparison
  template <class T>
  struct compare_astar
  {
    bool operator()(T a1, T a2) const
    {
      double f1 = a1->g + a1->h;
      double f2 = a2->g + a2->h;
      if( ( f1 >= f2 - 0.000001) && (f1 <= f2 +0.000001) )
        return a1->g < a2->g; // if equal compare gvals
      return f1 > f2;
    }
  };


  ///Define priority queue
  struct State; // forward declaration
  using StatePtr = std::shared_ptr<State>;
  using priorityQueue = boost::heap::d_ary_heap<StatePtr, boost::heap::mutable_<true>,
                        boost::heap::arity<2>, boost::heap::compare< compare_astar<StatePtr> >>;

  ///Lattice of the graph in graph search
  struct State
  {
    /// ID
    int id;
    /// Coord
    int x, y, z = 0;
    /// direction 
    int dx, dy, dz;                            // discrete coordinates of this node
    // id of predicessors
    int parentId = -1;

    // pointer to heap location
    priorityQueue::handle_type heapkey;

    // plan data
    double g = std::numeric_limits<double>::infinity();
    double h;
    bool opened = false;
    bool closed = false;

    /// 2D constructor
    State(int id, int x, int y, int dx, int dy )
      : id(id), x(x), y(y), dx(dx), dy(dy)
    {}

    /// 3D constructor
    State(int id, int x, int y, int z, int dx, int dy, int dz )
      : id(id), x(x), y(y), z(z), dx(dx), dy(dy), dz(dz)
    {}

  };


  /**
   * @brief GraphSearch class
   *
   * Implement A* and Jump Point Search
   */
  class GraphSearch
  {
    public:
     /**
       * @brief Graph search
       *
       */
      GraphSearch(const char* cMap, int xDim, int yDim,  double eps = 1, bool verbose = false);

      bool plan(int xStart, int yStart, int xGoal, int yGoal, int max_expand = -1);

      void useJps();
      std::vector<StatePtr> getPath() const;

      std::vector<StatePtr> getOpenedState() const;

    private:
      void getSucc(const StatePtr& curr, std::vector<int>& succ_ids, std::vector<double>& succ_costs);
      void getJpsSucc(const StatePtr& curr, std::vector<int>& succ_ids, std::vector<double>& succ_costs);
      std::vector<StatePtr> recoverPath(StatePtr node, int id);
      int coordToId(int x, int y) const;
      bool isFree(int x, int y) const;
      bool isOccupied(int x, int y) const;
      double get_heur(int x, int y) const;
      double get_heur(int x, int y, int z) const;

      bool hasForced(int x, int y, int dx, int dy);
      std::vector<std::array<int, 2>> prune(int x, int y, int dx, int dy);
      bool jump(int x, int y, int dx, int dy, int& new_x, int& new_y);

      const char* cMap_;
      priorityQueue pq_;
      std::vector<StatePtr> hm_;
      std::vector<bool> seen_;
      int xDim_, yDim_, zDim_;
      double eps_;
      bool verbose_;
      std::vector<StatePtr> path_;
      
      const char val_free_ = 0;
      int xGoal_, yGoal_;

      bool jps_ = false;
      std::vector<std::vector<std::array<int, 2>>> add_map_;
      std::vector<std::vector<std::array<int, 2>>> obs_map_;

 };
}
#endif
