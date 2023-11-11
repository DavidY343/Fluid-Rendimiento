//
// Created by manuel on 10/11/23.
//
#include <vector>
#include "gtest/gtest.h"
#include "sim/grid.hpp"


TEST(MiFuncionTest, Indice) {
  std::vector<double> vector = {0.3, 0.3};
  std::vector<particula> particulas = {};
  particulas.emplace_back(0, 0.5, 0.5, 0.5, 0, 0, 0, 0, 0, 0);
  grid malla = grid(vector, particulas);
  EXPECT_EQ(malla.obtener_indice(0, 0, 0), 0);
}
