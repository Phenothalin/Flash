#include "RR_vl.hpp"

#include "cluster_manager.hpp"
#include "component.hpp"
#include "data_warehouse_type.hpp"
#include "database_connection_pool_manager.hpp"
#include "material_object.hpp"
#include "model_parameter.hpp"

#include <cmath>
#include <gtest/gtest.h>
#include <memory>

using namespace material_object;
using namespace database;

class DatabaseTest : public ::testing::Test
{
protected:
  static void SetUpTestSuite() {
    if (!is_db_initialized)
    {
      std::cout << "Setting up database connection once for all tests...\n";
      initDatabase();
    }
  }

  static void TearDownTestSuite() {
    std::cout << "Cleaning up database connection...\n";
    // 这里可以添加清理数据库连接的代码（如果需要）
  }

  static auto isDatabaseInitialized() -> bool { return is_db_initialized; }

private:
  static bool is_db_initialized;
  static void initDatabase() {
    if (!DataManager::getInstance().initDataManager())
    {
      std::cout << "DataManager init failed!!!" << '\n';
      return;
    }

    if (!DataManager::getInstance().addDataWarehouses(
            DataWarehouseType::REMOTE_PUBLIC, true))
    {
      std::cout << "DatabaseConnectionPoolManager addDataWarehouses "
                   "REMOTE_PUBLIC failed!!!"
                << '\n';
      return;
    }

    if (!DataManager::getInstance().addDataWarehouses(
            DataWarehouseType::REMOTE_PRIVATE))
    {
      std::cout << "DatabaseConnectionPoolManager addDataWarehouses "
                   "REMOTE_PRIVATE failed!!!"
                << '\n';
      return;
    }

    is_db_initialized = true;
    std::cout << "Database initialized successfully.\n";
  }
};

class RRFlashTest : public ::DatabaseTest
{
protected:
  void TearDown() override {
    // 清理资源
    material_object::SubstanceManager::getInstance().clearAllSubstances();
  }
};
auto getManagerCluster(const std::vector<std::string> &substance_list)
    -> std::shared_ptr<material_object::Cluster> {
  // 1. 创建带有数据库加载物质的集群
  auto &substance_manager = material_object::SubstanceManager::getInstance();
  auto &cluster_manager = material_object::ClusterManager::getInstance();
  for (const auto &substance_name : substance_list)

  {
    substance_manager.createSubstanceFromDatabase(
        substance_name, database::DataWarehouseType::REMOTE_PRIVATE);
  }
  cluster_manager.createCluster("TestCluster");
  cluster_manager.addBatchSubstancesToCluster(
      "TestCluster", substance_manager.getAllSubstances());
  auto cluster = cluster_manager.getCluster("TestCluster");
  return cluster;
}
std::vector<string> feed1 = {"NITROGEN","CARBON DIOXIDE","METHANE",   "ETHANE",   "PROPANE",
                             "ISOBUTANE", "n-BUTANE", "ISOPENTANE",
                             "n-PENTANE", "n-HEXANE", "n-HEPTANE"};
std::vector<double> composition1 = {0.003,0.015,0.55, 0.14,  0.12,  0.05, 0.045,
                                    0.03,  0.025, 0.012, 0.01};
std::vector<string> feed2 = {"METHANE",   "ETHANE",   "PROPANE",
                              "ISOBUTANE", "n-BUTANE", "ISOPENTANE",
                              "n-PENTANE", "n-HEXANE", "n-HEPTANE","n-OCTANE"};
std::vector<double> composition2 = {0.55, 0.16,0.12,0.05,0.04,0.025,0.02,0.015,0.01,0.01};
TEST_F(RRFlashTest, PTFlash) {
    auto cluster = getManagerCluster(feed1);
    property_package::PropertyPackage SRKproperty(SRK, cluster);
    PTFlash ptflash(2e6, 295, composition1, SRKproperty);
    ptflash.calculate(ConvergenceMethod::NEWTON_RAPHSON);
    double temperature = ptflash.getTemperature();
    double pressure = ptflash.getPressure();
    double vaporfraction = ptflash.getVaporFraction();
    std::vector<double> molefraction_vapor = ptflash.getVapComp();
    std::vector<double> molefraction_liquid = ptflash.getLiqComp();
    std::cout << "T : " << temperature << " K\n"
              << "P : " << pressure << " Pa\n"
              << "vaporfraction " << vaporfraction << "\n";
    for (size_t i = 0; i < molefraction_vapor.size(); ++i)
    {
      std::cout << "组分" << i + 1 << "气相分率: " << molefraction_vapor[i]
                << "\n";
    }
    for (size_t i = 0; i < molefraction_liquid.size(); ++i)
    {
      std::cout << "组分" << i + 1 << "液相分率: " << molefraction_liquid[i]
                << "\n";
    }
}
bool DatabaseTest::is_db_initialized = false;

// TEST_F(RRFlashTest, PVFlash) {
//   auto cluster = getManagerCluster(feed1);
//   property_package::PropertyPackage SRKproperty(SRK, cluster);
//   PVFlash pvflash(5.5E6, 0.8, composition1, SRKproperty);
//   pvflash.calculate(ConvergenceMethod::HALLEY);
//   double temperature = pvflash.getTemperature();
//   double pressure = pvflash.getPressure();
//   double vaporfraction = pvflash.getVaporFraction();
//   std::vector<double> molefraction_vapor = pvflash.getVapComp();
//   std::vector<double> molefraction_liquid = pvflash.getLiqComp();
//   std::cout << "T : " << temperature << " K\n"
//             << "P : " << pressure << " Pa\n"
//             << "vaporfraction " << vaporfraction << "\n";
//   for (size_t i = 0; i < molefraction_vapor.size(); ++i)
//   {
//     std::cout << "组分" << i + 1 << "气相分率: " << molefraction_vapor[i]
//               << "\n";
//   }
//   for (size_t i = 0; i < molefraction_liquid.size(); ++i)
//   {
//     std::cout << "组分" << i + 1 << "液相分率: " << molefraction_liquid[i]
//               << "\n";
//   }
// }
// TEST_F(RRFlashTest, TVFlash) {
//   auto cluster = getManagerCluster(feed2);
//   property_package::PropertyPackage SRKproperty(SRK, cluster);
//   TVFlash tvflash(330, 0.7, composition2, SRKproperty);
//   tvflash.calculate(ConvergenceMethod::NEWTON_RAPHSON);
//   double temperature = tvflash.getTemperature();
//   double pressure = tvflash.getPressure();
//   double vaporfraction = tvflash.getVaporFraction();
//   std::vector<double> molefraction_vapor = tvflash.getVapComp();
//   std::vector<double> molefraction_liquid = tvflash.getLiqComp();
//   std::cout << "T : " << temperature << " K\n"
//             << "P : " << pressure << " Pa\n"
//             << "vaporfraction " << vaporfraction << "\n";
//   for (size_t i = 0; i < molefraction_vapor.size(); ++i)
//   {
//     std::cout << "组分" << i + 1 << "气相分率: " << molefraction_vapor[i]
//               << "\n";
//   }
//   for (size_t i = 0; i < molefraction_liquid.size(); ++i)
//   {
//     std::cout << "组分" << i + 1 << "液相分率: " << molefraction_liquid[i]
//               << "\n";
//   }
// }