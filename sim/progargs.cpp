//
// Created by david on 11/7/23.
//



//
// Created by david on 10/3/23.
//

#include "progargs.hpp"
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numbers>
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

// Función que comprueba los parámetros del programa, incluyendo el número de pasos de tiempo y la validez del formato
int comprobar_params(std::vector<std::string> const & argumentos, std::ifstream const & inputFile,
                     std::ofstream const & outputFile) {
  try {
    int const nSteps = std::stoi(argumentos[0]);
    if (nSteps < 0) {
      std::cerr << "Error: Invalid number of time steps."
                << "\n";
      return -2;
    }
  } catch (
      std::invalid_argument const & excepcion) {  // Solo ocurrira en caso de que argv[1] no sea int
    std::cerr << "Error: timesteps must be numeric."
              << "\n";
    return -1;
  }

  return comprobar_archivos(argumentos, inputFile, outputFile);
}


// Función que verifica el número de parámetros de la línea de comandos
int n_params(int argc) {
  if (argc != 4) {
    std::cerr << "Error: Invalid number of arguments:" << argc << "\n";
    return -1;
  }
  return 0;
}

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

//Comprobamos que los parametros que leemos de la cabecera estan bien
void comprobar_fallos_cabecera(std::vector<particula> const & particulas, int n_particulas_int) {
  int const longitud = static_cast<int>(particulas.size());
  int const error    = -5;
  if (n_particulas_int <= 0) {
    std::cerr << "Invalid number of particles: " << n_particulas_int << ".\n";
    std::exit(error);
  }
  if (longitud != n_particulas_int) {
    std::cerr << "Number of particles mismatch. Header: " << n_particulas_int
              << ", Found: " << longitud << ".\n";
    std::exit(error);
  }
}

// Con esta funcion creamos el array de partículas, y generamos los datos del archivo de entrada
std::vector<particula> crear_particulas(std::ifstream const & inputFile) {
  // Obtenemos el número de partículas de la cabecera
  auto n_particulas_int = read_binary_value<int>((std::istream &) inputFile);
  std::vector<particula> particulas;  // Creamos el vector con las partículas
  int ident = 0;
  while (!inputFile.eof()) {
    // Leer los datos de una partícula
    int const id_dat = ident;
    auto px_dat      = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto py_dat      = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto pz_dat      = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto hvx         = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto hvy         = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto hvz         = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto vx_dat      = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto vy_dat      = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    auto vz_dat      = static_cast<double>(read_binary_value<float>((std::istream &) inputFile));
    if (inputFile.eof()) { break; }
    particula const n_particula(id_dat, px_dat, py_dat, pz_dat, hvx, hvy, hvz, vx_dat, vy_dat,
                                vz_dat);
    particulas.push_back(n_particula);
    ident++;
  }
  comprobar_fallos_cabecera(particulas, n_particulas_int);
  return particulas;
}

// Inicializa la densidad y aceleración de la partícula
void particula::inicializar_densidad_aceleracion() {
  p  = 0;
  ax = constantes::g_const[0];
  ay = constantes::g_const[1];
  az = constantes::g_const[2];
}

// Interactúa con otra partícula para actualizar la densidad
void particula::interactuar_densidad(particula & part, double h) {
  double const distancia_cuadrado = calcularDistancia(part);

  if (distancia_cuadrado < (h * h)) {
    double const incremento_densidad = pow((h * h - distancia_cuadrado), 3);
    setp(getp() + incremento_densidad);
    double const incremento_densidad2 = pow((h * h - distancia_cuadrado), 3);
    part.setp(part.getp() + incremento_densidad2);
  }
}

// Calcula la distancia euclidiana al cuadrado entre dos partículas
double particula::calcularDistancia(particula const & part) const {
  double distancia = 0.0;

  // Calcula la distancia euclidiana entre las dos partículas
  double diff  = getpx() - part.getpx();
  distancia   += diff * diff;
  diff         = getpy() - part.getpy();
  distancia   += diff * diff;
  diff         = getpz() - part.getpz();
  distancia   += diff * diff;

  return distancia;  // En realidad es el modulo al cuadrado
}

// Interactúa con otra partícula para actualizar la aceleración
void particula::interactuar_aceleracion(particula & part, double h, double m) {
  double const distancia_cuadrado = calcularDistancia(part);
  std::vector<int> const r_num    = {15, 6, 45};
  if (distancia_cuadrado < pow(h, 2)) {
    double const dist       = sqrt(std::max(distancia_cuadrado, pow(10, -12)));
    std::vector<double> d_a = calcular_d_a(part, h, dist, m);
    interactuar_aceleracion2(part, d_a);
  }
}

// Calcula el factor de aceleración para la interacción entre partículas
double particula::calcular_factor_aceleracion(double h, double dist, double m,
                                              particula const & part) const {
  std::vector<int> const r_num = {15, 6};
  return (r_num[0] / (std::numbers::pi * pow(h, r_num[1]))) * ((3 * m * constantes::ps_const) / 2) *
         ((pow((h - dist), 2)) / dist) * (getp() + part.getp() - 2 * constantes::p_const);
}

// Calcula el factor de velocidad para la interacción entre partículas
double particula::calcular_factor_v(double h, double m) {
  std::vector<int> const r_num = {45, 6};
  return (r_num[0] / (std::numbers::pi * pow(h, r_num[1]))) * constantes::u_const * m;
}

// Calcula la aceleración resultante de la interacción entre partículas
std::vector<double> particula::calcular_d_a(particula const & part, double h, double dist,
                                            double m) const {
  double const factor_aceleracion = calcular_factor_aceleracion(h, dist, m, part);
  double const factor_v           = calcular_factor_v(h, m);

  std::vector<double> d_a = {
    (((getpx() - part.getpx()) * factor_aceleracion + (part.getvx() - getvx()) * factor_v) /
     (getp() * part.getp())),
    (((getpy() - part.getpy()) * factor_aceleracion + (part.getvy() - getvy()) * factor_v) /
     (getp() * part.getp())),
    (((getpz() - part.getpz()) * factor_aceleracion + (part.getvz() - getvz()) * factor_v) /
     (getp() * part.getp()))};
  return d_a;
}

// Actualiza la aceleración de la partícula basándose en la interacción con otra partícula
void particula::interactuar_aceleracion2(particula & part, std::vector<double> & d_a) {
  ax += d_a[0];
  ay += d_a[1];
  az += d_a[2];
  part.ax -= d_a[0];
  part.ay -= d_a[1];
  part.az -= d_a[2];
}

// Transforma la densidad de la partícula según una fórmula específica
void particula::transformar_densidad(double h, double m) {
  int const pow_6 = 6;
  int const d_315 = 315;
  int const d_64  = 64;
  int const d_9   = 9;
  p               = (p + pow(h, pow_6)) * (d_315 / (d_64 * std::numbers::pi * pow(h, d_9))) * m;
}

// Realiza la colisión con el límite en el eje X
void particula::colisionLimiteEjeX(bool lim_inf) {
  double const min_value = 0.0000000001;
  double const new_x     = getpx() + gethvx() * constantes::t_const;
  if (lim_inf) {
    double const difLimX = constantes::dp_const - (new_x - constantes::bmin_const[0]);
    if (difLimX > min_value) {
      setax(getax() + (constantes::sc_const * difLimX - constantes::dv_const * getvx()));
    }
  } else {
    double const difLimX = constantes::dp_const - (constantes::bmax_const[0] - new_x);
    if (difLimX > min_value) {
      setax(getax() - (constantes::sc_const * difLimX + constantes::dv_const * getvx()));
    }
  }
}

// Realiza la colisión con el límite en el eje Y
void particula::colisionLimiteEjeY(bool lim_inf) {
  double const min_value = 0.0000000001;
  double const new_y     = getpy() + gethvy() * constantes::t_const;
  if (lim_inf) {
    double const difLimY = constantes::dp_const - (new_y - constantes::bmin_const[1]);
    if (difLimY > min_value) {
      setay(getay() + (constantes::sc_const * difLimY - constantes::dv_const * getvy()));
    }
  } else {
    double const difLimY = constantes::dp_const - (constantes::bmax_const[1] - new_y);
    if (difLimY > min_value) {
      setay(getay() - (constantes::sc_const * difLimY + constantes::dv_const * getvy()));
    }
  }
}

// Realiza la colisión con el límite en el eje Z
void particula::colisionLimiteEjeZ(bool lim_inf) {
  double const min_value = 0.0000000001;
  double const new_z     = getpz() + gethvz() * constantes::t_const;
  if (lim_inf) {
    double const difLimZ = constantes::dp_const - (new_z - constantes::bmin_const[2]);
    if (difLimZ > min_value) {
      setaz(getaz() + (constantes::sc_const * difLimZ - constantes::dv_const * getvz()));
    }
  } else {
    double const difLimZ = constantes::dp_const - (constantes::bmax_const[2] - new_z);
    if (difLimZ > min_value) {
      setaz(getaz() - (constantes::sc_const * difLimZ + constantes::dv_const * getvz()));
    }
  }
}

// Realiza la colisión con el límite en la coordenada x del recinto
void particula::limiteRecintox(bool lim_inf) {
  double const dx_dat =
      lim_inf ? (getpx() - constantes::bmin_const[0]) : (constantes::bmax_const[0] - getpx());
  if (dx_dat < 0) {
    px  = lim_inf ? (constantes::bmin_const[0] - dx_dat) : (constantes::bmax_const[0] + dx_dat);
    vx  = -vx;
    hvx = -hvx;
  }
}

// Realiza la colisión con el límite en la coordenada Y del recinto
void particula::limiteRecintoy(bool lim_inf) {
  double const dy_dat =
      lim_inf ? (getpy() - constantes::bmin_const[1]) : (constantes::bmax_const[1] - getpy());
  if (dy_dat < 0) {
    py  = lim_inf ? (constantes::bmin_const[1] - dy_dat) : (constantes::bmax_const[1] + dy_dat);
    vy  = -vy;
    hvy = -hvy;
  }
}

// Realiza la colisión con el límite en la coordenada Z del recinto
void particula::limiteRecintoz(bool lim_inf) {
  double const dz_dat =
      lim_inf ? (getpz() - constantes::bmin_const[2]) : (constantes::bmax_const[2] - getpz());
  if (dz_dat < 0) {
    pz  = lim_inf ? (constantes::bmin_const[2] - dz_dat) : (constantes::bmax_const[2] + dz_dat);
    vz  = -vz;
    hvz = -hvz;
  }
}

// Actualiza el movimiento de la partícula según la simulación
void particula::actualizarMovimiento() {
  setpx(getpx() + gethvx() * constantes::t_const + getax() * pow(constantes::t_const, 2));
  setvx(gethvx() + (getax() * constantes::t_const) / 2);
  sethvx(gethvx() + getax() * constantes::t_const);

  setpy(getpy() + gethvy() * constantes::t_const + getay() * pow(constantes::t_const, 2));
  setvy(gethvy() + (getay() * constantes::t_const) / 2);
  sethvy(gethvy() + getay() * constantes::t_const);

  setpz(getpz() + gethvz() * constantes::t_const + getaz() * pow(constantes::t_const, 2));
  setvz(gethvz() + (getaz() * constantes::t_const) / 2);
  sethvz(gethvz() + getaz() * constantes::t_const);
}