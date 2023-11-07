//
// Created by david on 11/7/23.
//

#include "progargs.hpp"
//
// Created by david on 10/3/23.
//

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// Funcion que comprueba que los archivos se abrieron correctamente
int comprobar_archivos(std::vector<std::string> const & argumentos, std::ifstream const & inputFile,
                       std::ofstream const & outputFile) {
  if (!inputFile.is_open()) {  // Si nos devuelve False, error
    std::cerr << "Error: Cannot open " << argumentos[1] << " for reading" << '\n';
    return -3;
  }
  if (!outputFile.is_open()) {  // Si nos devuelve False, error
    std::cerr << "Error: Cannot open " << argumentos[2] << " for writing" << '\n';
    return -4;
  }
  return 0;
}

// Funcion que comprueba el número de parámetros pasados y si iteraciones es un numero
int comprobar_params(int argc, std::vector<std::string> const & argumentos,
                     std::ifstream const & inputFile, std::ofstream const & outputFile) {
  if (argc != 4) {
    std::cerr << "Error: Invalid number of arguments:" << argc << "\n";
    return -1;
  }
  int nSteps = 0;
  try {
    nSteps = std::stoi(argumentos[0]);
  } catch (
      std::invalid_argument const & excepcion) {  // Solo ocurrira en caso de que argv[1] no sea int
    std::cerr << "Error: timesteps must be numeric."
              << "\n";
    return -1;
  }
  if (nSteps < 0) {
    std::cerr << "Error: Invalid number of timesteps."
              << "\n";
    return -2;
  }
  return comprobar_archivos(argumentos, inputFile, outputFile);
};

// Funciones para leer del archivo
template <typename T>
  requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
char * as_writable_buffer(T & value) {
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  return reinterpret_cast<char *>(&value);
}

template <typename T>
  requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
char const * as_buffer(T const & value) {
  // NOLINTNEXTLINE(cppcoreguidelines-pro-type-reinterpret-cast)
  return reinterpret_cast<char const *>(&value);
}

template <typename T>
  requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
T read_binary_value(std::istream & is) {
  T value{};
  is.read(as_writable_buffer(value), sizeof(value));
  return value;
}

std::vector<double> calculo_tmñ_malla(std::vector<double> const & max_min, double h) {
  double const nx_vect = (max_min[0] - max_min[1]) / h;
  double const ny_vect = (max_min[2] - max_min[3]) / h;
  double const nz_vect = (max_min[4] - max_min[5]) / h;

  return {nx_vect, ny_vect, nz_vect};
}

std::vector<double> calcular_max_min(std::vector<particula> const & particulas) {
  double max_px = -std::numeric_limits<double>::max();
  double min_px = std::numeric_limits<double>::max();
  double max_py = -std::numeric_limits<double>::max();
  double min_py = std::numeric_limits<double>::max();
  double max_pz = -std::numeric_limits<double>::max();
  double min_pz = std::numeric_limits<double>::max();

  for (particula const & n_particula : particulas) {
    double const px_dat = n_particula.getpx();
    double const py_dat = n_particula.getpy();
    double const pz_dat = n_particula.getpz();

    max_px = std::max(max_px, px_dat);
    min_px = std::min(min_px, px_dat);
    max_py = std::max(max_py, py_dat);
    min_py = std::min(min_py, py_dat);
    max_pz = std::max(max_pz, pz_dat);
    min_pz = std::min(min_pz, pz_dat);
  }

  std::vector const max_min{max_px, min_px, max_py, min_py, max_pz, min_pz};
  return max_min;
}

std::vector<double> calculo_tmñ_bloque_malla(std::vector<double> const & max_min,
                                             std::vector<double> const & n_const) {
  double const sx_vect = (max_min[0] - max_min[1]) / n_const[0];
  double const sy_vect = (max_min[2] - max_min[3]) / n_const[1];
  double const sz_vect = (max_min[4] - max_min[5]) / n_const[2];

  return {sx_vect, sy_vect, sz_vect};
}

std::vector<double> longitud_masa(std::ifstream & inputFile) {
  inputFile.seekg(0, std::ios::beg);  // situamos el puntero del archivo en la cabecera
  auto ppm = static_cast<double>(read_binary_value<float>(inputFile));

  // Realizamos los cálculos necesarios
  double const m_dat = constantes::p_const / std::pow(ppm, 3.0);
  double const h_dat = constantes::r_const / ppm;

  return {m_dat, h_dat};
}

int comprobar_fallos_cabecera(std::vector<particula> const & particulas,
                              std::ifstream & inputFile) {
  inputFile.seekg(1, std::ios::beg);  // situamos el puntero del archivo
  // Obtenemos el número de partículas de la cabecera
  auto n_particulas_int = read_binary_value<int>(inputFile);
  size_t const longitud = particulas.size();
  if (longitud != 0) {
    int const error = -5;
    if (longitud <= 0) {
      std::cerr << "Invalid number of particles:" << longitud;
    } else {
      std::cerr << "Number of particles mismatch. Header:" << n_particulas_int
                << "Found:" << longitud;
    }
    return error;
  }
  return 0;
}

// Con esta funcion creamos el array de partículas, y generamos los datos del archivo de entrada
std::pair<int, std::vector<particula>> crear_particulas(std::ifstream & inputFile) {
  inputFile.seekg(2, std::ios::beg);  // nos aseguramos de empezar después de la cabecera
  std::vector<particula> particulas;  // Creamos el vector con las partículas
  int ident = 0;

  while (!inputFile.eof()) {
    // Leer los datos de una partícula
    int const id_dat = ident;
    auto px_dat      = static_cast<double>(read_binary_value<float>(inputFile));
    auto py_dat      = static_cast<double>(read_binary_value<float>(inputFile));
    auto pz_dat      = static_cast<double>(read_binary_value<float>(inputFile));
    auto hvx         = static_cast<double>(read_binary_value<float>(inputFile));
    auto hvy         = static_cast<double>(read_binary_value<float>(inputFile));
    auto hvz         = static_cast<double>(read_binary_value<float>(inputFile));
    auto vx_dat      = static_cast<double>(read_binary_value<float>(inputFile));
    auto vy_dat      = static_cast<double>(read_binary_value<float>(inputFile));
    auto vz_dat      = static_cast<double>(read_binary_value<float>(inputFile));
    particula const n_particula(id_dat, px_dat, py_dat, pz_dat, hvx, hvy, hvz, vy_dat, vx_dat,
                                vz_dat);
    particulas.push_back(n_particula);
    ident++;
  }
  int const error = comprobar_fallos_cabecera(particulas, inputFile);
  if (error < 0) { return std::make_pair(error, std::vector<particula>()); }
  return std::make_pair(0, particulas);
}

void  init_params(std::ifstream & inputFile)
{
  const std::vector<double> cabeceras = longitud_masa(inputFile);

}