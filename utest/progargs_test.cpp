//
// Created by albavidales on 11/13/23.
//

#include "sim/progargs.hpp"

#include "gtest/gtest.h"
#include <fstream>

TEST(NParamsTest, NumeroArgumentosValido) {
  // Caso de prueba: Número de argumentos válido
  EXPECT_EQ(0, n_params(4));
}

TEST(NParamsTest, NumeroArgumentosInvalido) {
  // Caso de prueba: Número de argumentos inválido
  EXPECT_EQ(-1, n_params(3));
}

TEST(ComprobarParamsTest, ParamsCorrectos) {
  // Caso de prueba: todos los argumentos son correctos
  std::vector<std::string> const args{"10", "small.fld", "output.fld"};
  std::ifstream const inputFile(args[1],
                                std::ios::binary);  // Creamos un objeto inputFile para lecutura
  std::ofstream const outputFile(args[2], std::ios::binary);
  EXPECT_EQ(0, comprobar_params(args, inputFile, outputFile));
}

TEST(ComprobarParamsTest, TimeStepsNoNumerico) {
  // Caso de prueba: Timesteps no numéricos
  std::vector<std::string> const args{"ey", "small.fld", "output.fld"};
  std::ifstream const inputFile(args[1],
                                std::ios::binary);  // Creamos un objeto inputFile para lecutura
  std::ofstream const outputFile(args[2], std::ios::binary);
  EXPECT_EQ(-1, comprobar_params(args, inputFile, outputFile));
}

TEST(ComprobarParamsTest, TimeStepsNegativos) {
  // Caso de prueba: Timesteps negativos
  std::vector<std::string> const args{"-5", "small.fld", "output.fld"};
  std::ifstream const inputFile(args[1],
                                std::ios::binary);  // Creamos un objeto inputFile para lecutura
  std::ofstream const outputFile(args[2], std::ios::binary);
  EXPECT_EQ(-2, comprobar_params(args, inputFile, outputFile));
}

TEST(ComprobarArchivosTest, ErrorAlAbrirInputFile) {
  // Caso de prueba: Input no se abre correctamente
  std::vector<std::string> const args{"10", "no_existe_input.fld", "output.fld"};
  std::ifstream const inputFile(args[1],
                                std::ios::binary);  // Creamos un objeto inputFile para lecutura
  std::ofstream const outputFile(args[2], std::ios::binary);
  EXPECT_EQ(-3, comprobar_params(args, inputFile, outputFile));
}

// PROBLEMA PORQUE SIEMPRE ME LO CREA XD
TEST(ComprobarArchivosTest, ErrorAlAbrirOutputFile) {
  // Caso de prueba: Output no se abre correctamente
  std::vector<std::string> const args{"10", "small.fld", "no_existe_output.txt"};
  std::ifstream const inputFile(args[1],
                                std::ios::binary);  // Creamos un objeto inputFile para lecutura
  std::ofstream const outputFile(args[2], std::ios::binary);
  EXPECT_EQ(-4, comprobar_params(args, inputFile, outputFile));
}