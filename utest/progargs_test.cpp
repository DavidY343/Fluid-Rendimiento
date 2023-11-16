//
// Created by albavidales on 11/13/23.
//

#include "sim/progargs.hpp"

#include "gtest/gtest.h"
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

// EN ESTOS TRES HABÍA QUE CAMBIAR LO DEL ppm
TEST(ComprobarFallosCabeceraTest, InvalidNumberOfParticles) {
  // Caso de prueba: Número de partículas no válido
  std::ifstream const inputFile("small.fld",
                                std::ios::binary);  // Creamos un objeto inputFile para lectura
  auto ppm = static_cast<double>(
      read_binary_value<float>((std::istream &) inputFile));

  std::cout << "ppm es " << ppm << "\n";

  std::vector<particula> const particulas = crear_particulas(inputFile);

  int const n_particulas_int = 0;  // Número de partículas no válido

  ASSERT_EXIT(
      {
        // Llamamos a la función con un valor de partículas no válido
        comprobar_fallos_cabecera(particulas, n_particulas_int);
      },
      ::testing::ExitedWithCode(256-5), // Código de salida esperado
      "Invalid number of particles: 0."
  );
}

TEST(ComprobarFallosCabeceraDeathTest2, InvalidNumberOfParticles) {
  // Caso de prueba: Número de partículas diferente al de la cabecera
  std::ifstream const inputFile("small.fld",
                                std::ios::binary);  // Creamos un objeto inputFile para lectura
  auto ppm = static_cast<double>(
      read_binary_value<float>((std::istream &) inputFile));

  std::cout << "ppm es " << ppm << "\n";

  std::vector<particula> const particulas = crear_particulas(inputFile);

  int const n_particulas_int = 4799;  // Número de partículas no válido

  ASSERT_EXIT(
      {
        // Llamamos a la función con un valor de partículas no válido
        comprobar_fallos_cabecera(particulas, n_particulas_int);
      },
      ::testing::ExitedWithCode(256-5), // Código de salida esperado
      "Number of particles mismatch. Header: 4799, Found: 4800."
  );
}

TEST(ComprobarFallosCabeceraCorrecto, NumeroValido){
  // Caso de prueba: Numero de particulas creadas igual al de la cabecera
  std::ifstream const inputFile("small.fld",
                                std::ios::binary);  // Creamos un objeto inputFile para lectura
  auto ppm = static_cast<double>(
      read_binary_value<float>((std::istream &) inputFile));

  std::cout << "ppm es " << ppm << "\n";

  std::vector<particula> const particulas = crear_particulas(inputFile);

  int const n_particulas_int = 4800;  // Número de partículas correcto

  // Capturar la salida estándar para verificar el mensaje de error
  testing::internal::CaptureStderr();

  // Llamamos a la función con un valor de partículas válido
  comprobar_fallos_cabecera(particulas, n_particulas_int);

  // Verificar que la salida estándar esté vacía (sin mensaje de error)
  std::string const output = testing::internal::GetCapturedStderr();
  EXPECT_TRUE(output.empty());
}
