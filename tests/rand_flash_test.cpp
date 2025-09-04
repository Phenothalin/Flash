// rand_flash/tests/rand_flash_test.cpp
#include "rand_flash.hpp"

#include "cluster_manager.hpp"
#include "component.hpp"
#include "data_warehouse_type.hpp"
#include "database_connection_pool_manager.hpp"
#include "material_object.hpp"
#include "model_parameter.hpp"

#include <gtest/gtest.h>
#include <memory>
#include <cmath>
constexpr double MOLAR_WEIGHT_PROPANE          = 44.097;     // g/mol
constexpr double CRITICAL_TEMPERATURE_PROPANE  = 369.83;     // K
constexpr double ACENTRIC_FACTOR_PROPANE       = 0.152;      // 无量纲
constexpr double CRITICAL_PRESSURE_PROPANE     = 4.248e6;    // Pa (≈42.48 bar)

// Test for creating PropertyPackage with IDEAL model
constexpr double TEMPERATURE = 295.024;
constexpr double PRESSURE = 1832850;
constexpr double MOLE_FRACTION_WATER_LIQUID = 0.532813581049914;
constexpr double MOLE_FRACTION_ETHANE_LIQUID = 0.467186418950086;
constexpr double MOLE_FRACTION_WATER_VAPOR = 0.000763158738719824;
constexpr double MOLE_FRACTION_ETHANE_VAPOR = 0.99923684126128;
//  water
constexpr double CRITICAL_TEMPERATURE_WATER = 647.096;
constexpr double CRITICAL_PRESSURE_WATER = 22064000;
constexpr double ACENTRIC_FACTOR_WATER = 0.229;
constexpr double MOLAR_WEIGHT_WATER = 18.01528;
// ethane
constexpr double CRITIAL_TEMPERATURE_ETHANE = 305.32;


constexpr double CRITIAL_PRESSURE_ETHANE = 4872000;
constexpr double ACENTRIC_FACTOR_ETHANE = 0.099493;
constexpr double MOLAR_WEIGHT_ETHANE = 30.07;

using namespace randflash;
using namespace thermo;
using namespace material_object;
using namespace database;

auto getTestCluster() -> std::shared_ptr<material_object::Cluster> {
  material_object::Cluster cluster("cluster");
  //
  // —— 第一组：乙烷（Ethane）
  //
  const std::vector<double> PARA_ETHANE{
      51.857,     // DIPPR101 a
      -2598.7,    // DIPPR101 b
      -5.1283,    // DIPPR101 c
      0.000014913,// DIPPR101 d
      2,          // DIPPR101 e
      100,        // Tmin
      500         // Tmax
  };
  const std::vector<double> MOLAR_DENSITY_PARA_ETHANE{
      1.9122,  // DIPPR105 A
      0.27937, // DIPPR105 B
      305.32,  // DIPPR105 C
      0.29187  // DIPPR105 D
  };
  material_object::Substance ethane("Ethane", "74-84-0");
  // 常数属性
  ethane.setConstantProperty(database::ConstantPropertyType::MOLECULAR_WEIGHT,
                             MOLAR_WEIGHT_ETHANE);
  ethane.setConstantProperty(database::ConstantPropertyType::CRITICAL_TEMPERATURE,
                             CRITIAL_TEMPERATURE_ETHANE);
  ethane.setConstantProperty(database::ConstantPropertyType::ACENTRIC_FACTOR,
                             ACENTRIC_FACTOR_ETHANE);
  ethane.setConstantProperty(database::ConstantPropertyType::CRITICAL_PRESSURE,
                             CRITIAL_PRESSURE_ETHANE);
  // 温度依赖属性
  ethane.setTemperatureDependentProperty(
      database::TemperaturePropertyType::VAPOR_PRESSURE,
      database::EquationType::DIPPR101,
      PARA_ETHANE,
      100.0, 500.0);
  ethane.setTemperatureDependentProperty(
      database::TemperaturePropertyType::DENSITY_OF_LIQUID,
      database::EquationType::DIPPR105,
      MOLAR_DENSITY_PARA_ETHANE,
      100.0, 500.0);
  ethane.setTemperatureDependentProperty(
      database::TemperaturePropertyType::IDEAL_GAS_HEAT_CAPACITY,
      database::EquationType::DIPPR127,
      std::vector<double>{33257.8886, 70253.576162, -1605.767453,
                          60622.500697, -3784.682119, 12107.569262,
                          489.57675},
      20.0, 1500.0);
  ethane.setTemperatureDependentProperty(
      database::TemperaturePropertyType::HEAT_CAPACITY_OF_LIQUID,
      database::EquationType::DIPPR114,
      std::vector<double>{44.009, 89718.0, 918.77, -1886.0},
      92.0, 290.0,
      CRITIAL_TEMPERATURE_ETHANE);
  ethane.setTemperatureDependentProperty(
      database::TemperaturePropertyType::HEAT_OF_VAPORIZATION,
      database::EquationType::DIPPR106,
      std::vector<double>{21091000.0, 0.60646, -0.55492, 0.32799},
      90.35, 305.32,
      CRITIAL_TEMPERATURE_ETHANE);
  // 把乙烷加到 Cluster
  material_object::Component component;
  component.setSubstance(std::make_shared<material_object::Substance>(ethane));
  cluster.addComponent(std::make_shared<material_object::Component>(component));
  //
  // —— 第二组：丙烷（Propane）
  //
  const std::vector<double> PARA_PROPANE{
      75.05,      // DIPPR101 a (示例值，请替换为准确文献/数据库值)
      -2200.3,    // DIPPR101 b
      -6.788,     // DIPPR101 c
      0.0000168,  // DIPPR101 d
      2,          // DIPPR101 e
      85.5,       // Tmin（丙烷凝固点约85.5 K）
      369.83      // Tmax（临界温度）
  };
  const std::vector<double> MOLAR_DENSITY_PARA_PROPANE{
      2.0098, // DIPPR105 A（示例）
      0.240,  // DIPPR105 B
      231.1,  // DIPPR105 C
      0.308   // DIPPR105 D
  };
  material_object::Substance propane("Propane", "74-98-6");
  // 常数属性
  propane.setConstantProperty(database::ConstantPropertyType::MOLECULAR_WEIGHT,
                              MOLAR_WEIGHT_PROPANE);
  propane.setConstantProperty(database::ConstantPropertyType::CRITICAL_TEMPERATURE,
                              CRITICAL_TEMPERATURE_PROPANE);
  propane.setConstantProperty(database::ConstantPropertyType::ACENTRIC_FACTOR,
                              ACENTRIC_FACTOR_PROPANE);
  propane.setConstantProperty(database::ConstantPropertyType::CRITICAL_PRESSURE,
                              CRITICAL_PRESSURE_PROPANE);
  // 温度依赖属性
  propane.setTemperatureDependentProperty(
      database::TemperaturePropertyType::VAPOR_PRESSURE,
      database::EquationType::DIPPR101,
      PARA_PROPANE,
      85.5, 369.83);
  propane.setTemperatureDependentProperty(
      database::TemperaturePropertyType::DENSITY_OF_LIQUID,
      database::EquationType::DIPPR105,
      MOLAR_DENSITY_PARA_PROPANE,
      85.5, 369.83);
  propane.setTemperatureDependentProperty(
      database::TemperaturePropertyType::IDEAL_GAS_HEAT_CAPACITY,
      database::EquationType::DIPPR127,
      std::vector<double>{34500.0, 70000.0, -1500.0,
                          60000.0, -4000.0, 12000.0,
                          500.0},   // 示例
      20.0, 1500.0);
  propane.setTemperatureDependentProperty(
      database::TemperaturePropertyType::HEAT_CAPACITY_OF_LIQUID,
      database::EquationType::DIPPR114,
      std::vector<double>{35.712, 105000.0, 900.0, -1200.0},  // 示例
      230.9, 369.83,
      CRITICAL_TEMPERATURE_PROPANE);
  propane.setTemperatureDependentProperty(
      database::TemperaturePropertyType::HEAT_OF_VAPORIZATION,
      database::EquationType::DIPPR106,
      std::vector<double>{22000000.0, 0.65, -0.6, 0.35},  // 示例
      230.9, 369.83,
      CRITICAL_TEMPERATURE_PROPANE);
  // 把丙烷加到 Cluster
  component.setSubstance(std::make_shared<material_object::Substance>(propane));
  cluster.addComponent(std::make_shared<material_object::Component>(component));

  return std::make_shared<material_object::Cluster>(cluster);
}

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


class RandFlashTest : public ::DatabaseTest {
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

// TEST_F(RandFlashTest, two_phase_database1){
//   auto cluster = getManagerCluster({"ETHANE", "PROPANE"});
//   randflash::RandFlash randflashSolver(
//       PropertyPackageType::SRK, cluster,
//       *ls::createEigenSolver());
//   double P = 1832850, T = 295;
//   std::vector<double> feed = {1.0, 1.0};
//   std::vector<std::vector<double>> elementMatrix(2, std::vector<double>(2, 0));
//   // C 原子
//   elementMatrix[0][0] = 2;
//   elementMatrix[0][1] = 3; 
//   // H 原子
//   elementMatrix[1][0] = 6;
//   elementMatrix[1][1] = 8;

//   FlashResult res = randflashSolver.solveTwoPhase(
//       P, T, feed, elementMatrix, {0.1, 0.9},{}, 10);
//   std::cout << "Convergence error: " << res.convergenceError << std::endl;
//   ASSERT_TRUE(res.success) << "Expected convergence";
//   std::cout << "vaporFraction: " << res.vaporFraction << std::endl;
//   std::cout << "vaporComposition: ";
//   for (const auto &comp : res.vaporComposition) {
//     std::cout << comp << " ";
//   }
//   std::cout << std::endl;
//   std::cout << "liquidComposition: ";
//   for (const auto &comp : res.liquidComposition) {
//     std::cout << comp << " "; 
//   }
//   std::cout << std::endl;
// }

TEST_F(RandFlashTest, two_phase_database2){
  auto cluster = getManagerCluster({"ETHANE", "ETHYLENE","PROPANE","PROPYLENE"});
  randflash::RandFlash randflashSolver(
      PropertyPackageType::SRK, cluster,
      *ls::createEigenSolver());
  double P = 2e6, T = 295;
  std::vector<double> feed = {1.0, 1.0, 1.0, 1.0};
  std::vector<std::vector<double>> elementMatrix(4, std::vector<double>(4, 0));
  // C 原子
  elementMatrix[0][0] = 1;
  elementMatrix[0][1] = 0; 
  elementMatrix[0][2] = 0; 
  elementMatrix[0][3] = 0; 
  // H 原子
  elementMatrix[1][0] = 0;
  elementMatrix[1][1] = 1;
  elementMatrix[1][2] = 0;
  elementMatrix[1][3] = 0; 

  elementMatrix[2][0] = 0;
  elementMatrix[2][1] = 0;
  elementMatrix[2][2] = 1;
  elementMatrix[2][3] = 0; 

  elementMatrix[3][0] = 0;
  elementMatrix[3][1] = 0;
  elementMatrix[3][2] = 0;
  elementMatrix[3][3] = 1; 
  FlashResult res = randflashSolver.solveTwoPhase(
      P, T, feed, elementMatrix, {0.2,0.3,0.2,0.3},{}, 10);
  std::cout << "Convergence error: " << res.convergenceError << std::endl;
  ASSERT_TRUE(res.success) << "Expected convergence";
  std::cout << "vaporFraction: " << res.vaporFraction << std::endl;
  std::cout << "vaporComposition: ";
  for (const auto &comp : res.vaporComposition) {
    std::cout << comp << " ";
  }
  std::cout << std::endl;
  std::cout << "liquidComposition: ";
  for (const auto &comp : res.liquidComposition) {
    std::cout << comp << " "; 
  }
  std::cout << std::endl;
}

bool DatabaseTest::is_db_initialized = false;