#include "sim/block.hpp"
#include "sim/grid.hpp"

#include "gtest/gtest.h"

TEST(DevolverParticulasTest, ComprobarDevolverParticulasBase) {
  double const m                    = 0.00011779;
  double const h                    = 0.00830882;
  std::vector<double> const vector  = {m, h};
  std::vector<particula> particulas = {};

  particula const part = particula(0, -0.1, -0.1, -0.1, 0, 0, 0, 0, 0, 0);
  particulas.push_back(part);
  grid malla = grid(vector, particulas);

  std::vector<block> bloques = malla.getbloques();

  EXPECT_EQ(bloques[0].devolver_particulas().size(), 1);
}