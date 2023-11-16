//
// Created by albavidales on 11/13/23.
//

#include "sim/progargs.hpp"

#include "gtest/gtest.h"
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>

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

TEST(ComprobarArchivosTest, ErrorAlAbrirOutputFile) {
  // Caso de prueba: Output no se abre correctamente
  std::vector<std::string> const args{"10", "small.fld"};
  std::ifstream const inputFile(args[1],
                                std::ios::binary);  // Creamos un objeto inputFile para lecutura
  std::ofstream const outputFile(args[2], std::ios::binary);  // Da error porque no hay nada
  EXPECT_EQ(-4, comprobar_params(args, inputFile, outputFile));
}

// PROBLEMA CON ESTE TEST, PUEDE  QUE HAYA QUE CAMBIAR EL EXIT, NO SER PUEDEN HACER TESTS PORQUE
// TERMINA CON EL FLUJO DEL PROGRAMA
TEST(ComprobarFallosCabeceraTest, InvalidNumberOfParticles) {
  // Caso de prueba: Número de partículas no válido
  std::ifstream const inputFile("small.fld",
                                std::ios::binary);  // Creamos un objeto inputFile para lectura
  std::vector<particula> const particulas = crear_particulas(inputFile);  // Agregamos punto y coma

  int const n_particulas_int = -2;  // Número de partículas no válido

  // Capturamos la salida estándar para verificar el mensaje de error
  testing::internal::CaptureStderr();

  // Llamamos a la función con un valor de partículas no válido
  comprobar_fallos_cabecera(particulas, n_particulas_int);

  // Verificamos el mensaje de error
  std::string const output = testing::internal::GetCapturedStderr();
  EXPECT_EQ(output, "Invalid number of particles: -2.\n");
}
