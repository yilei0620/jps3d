/**
 * @file astar_2d_util.h
 * @brief JPS::AStar2DUtil Class
 */

#include <jps3d/planner_base.h>
#include <jps3d/planner/graph_search.h>

namespace JPS {
  class AStar2DUtil: public PlannerBase {
    public:
      AStar2DUtil(bool verbose = false);
      bool plan(const Vec3f &start, const Vec3f &goal, decimal_t eps = 1.0);
    private:
      Vec3i _start_int;
      Vec3i _goal_int;
  };
}
