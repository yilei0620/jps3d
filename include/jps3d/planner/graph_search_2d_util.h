/**
 * @file graph_search_2d_util.h
 * @brief JPS::GraphSearch2DUtil Class
 */

#include <jps3d/planner_base.h>
#include <jps3d/planner/graph_search.h>

namespace JPS {
  class GraphSearch2DUtil: public PlannerBase {
    public:
      GraphSearch2DUtil(bool verbose = false);
      bool plan(const Vec3f &start, const Vec3f &goal, decimal_t eps = 1.0, bool use_jps = false);
  };
}
