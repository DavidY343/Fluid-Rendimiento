//
// Created by david on 11/7/23.
//
//
// Created by david on 10/3/23.
//

#include "grid.hpp"

#include "progargs.hpp"
#include "progargs.cpp"
#include <iostream>
#include <tuple>
/*
// Función para colisiones con límites en el eje x
void interaccionesLimitesEjeX(Particula& particula, int cx, double xmin, double xmax) {
  if (cx == 0 || cx == nx - 1) {
    double dx = (cx == 0) ? particula.px - xmin : xmax - particula.px;
    if (dx < 0) {
      particula.px = (cx == 0) ? xmin - dx : xmax + dx;
      particula.vx = -particula.vx;
      particula.hvx = -particula.hvx;
    }
  }
}

// Función para colisiones con límites en el eje y
void interaccionesLimitesEjeY(Particula& particula, int cy, double ymin, double ymax) {
  if (cy == 0 || cy == ny - 1) {
    double dy = (cy == 0) ? particula.py - ymin : ymax - particula.py;
    if (dy < 0) {
      particula.py = (cy == 0) ? ymin - dy : ymax + dy;
      particula.vy = -particula.vy;
      particula.hvy = -particula.hvy;
    }
  }
}

// Función para colisiones con límites en el eje z
void interaccionesLimitesEjeZ(Particula& particula, int cz, double zmin, double zmax) {
  if (cz == 0 || cz == nz - 1) {
    double dz = (cz == 0) ? particula.pz - zmin : zmax - particula.pz;
    if (dz < 0) {
      particula.pz = (cz == 0) ? zmin - dz : zmax + dz;
      particula.vz = -particula.vz;
      particula.hvz = -particula.hvz;
    }
  }
}

// Función general que llama a las tres funciones anteriores
void interaccionesLimitesRecinto(Particula& particula, int cx, int cy, int cz, double xmin, double
xmax, double ymin, double ymax, double zmin, double zmax) { interaccionesLimitesEjeX(particula, cx,
xmin, xmax); interaccionesLimitesEjeY(particula, cy, ymin, ymax);
  interaccionesLimitesEjeZ(particula, cz, zmin, zmax);
}

*/
void grid::colisiones_particulas() {
  std::vector<int> coordenadas;
  for (int i = 0; i < getnx() * getny() * getnz(); i++) {
    coordenadas = obtener_coordenadas(i);
    if (coordenadas[0] == 0 || coordenadas[0] == getnx() - 1) {
      bucle_colisiones(i, coordenadas[0] == 0, 0);
    } else if (coordenadas[1] == 0 || coordenadas[1] == getny() - 1) {
      bucle_colisiones(i, coordenadas[1] == 0, 1);
    } else if (coordenadas[2] == 0 || coordenadas[2] == getnz() - 1) {
      bucle_colisiones(i, coordenadas[1] == 0, 2);
    }
  }
}

void grid::bucle_colisiones(int num_bloque, bool lim_inf, int dimension) {
  for (auto & particula : bloques[num_bloque].particulas) {
      if (dimension == 0) {
        particula.colisionLimiteEjeX(lim_inf);
        // 4.3.4
        particula.actualizarMovimiento();
        // 4.3.5
        particula.limiteRecintox(lim_inf);
      }else if (dimension == 1) {
        particula.colisionLimiteEjeY(lim_inf);
        // 4.3.4
        particula.actualizarMovimiento();
        // 4.3.5
        particula.limiteRecintoy(lim_inf);
      }else if (dimension == 2) {
        particula.colisionLimiteEjeZ(lim_inf);
        // 4.3.4
        particula.actualizarMovimiento();
        // 4.3.5
        particula.limiteRecintoz(lim_inf);
      }
  }
}

grid init_params(std::ifstream const & inputFile) {
  std::vector<double> vtor                = longitud_masa(inputFile);
  std::vector<particula> const particulas = crear_particulas(inputFile);
  using namespace std;
  cout << "Inizializando malla con m=" << vtor[0] << " y h=" << vtor[1] << "\n";
  grid const malla(vtor, particulas);
  return malla;
}

void init_simulate(int const max_iteraciones, grid & malla) {
  using namespace std;
  for (int iteracion = 1; iteracion <= max_iteraciones; iteracion++) {
    cout << "****************************************************" << endl;
    cout << "iniciando  iteracion " << iteracion << endl;

    malla.simular();

    cout << "finalizada iteracion " << iteracion << endl;
    cout << "----------------------------------------------------" << endl;
  }
}
void escribir_parametros_generales(std::ofstream & outputFile, std::ifstream const & inputFile) {
  auto ppm = (read_binary_value<float>((std::istream &) inputFile));
  auto n_particulas_int = read_binary_value<int>((std::istream &) inputFile);
  outputFile.write(as_buffer(ppm), sizeof(ppm));
  outputFile.write(as_buffer(n_particulas_int), sizeof(n_particulas_int));
}

// Función para convertir los datos de las partículas de doble a simple precisión
std::tuple<float, float, float, float, float, float, float, float, float>
    convertirDatos(particula const & particula) {
  auto px_dat = static_cast<float>(particula.getpx());
  auto py_dat = static_cast<float>(particula.getpy());
  auto pz_dat = static_cast<float>(particula.getpz());
  auto hvx    = static_cast<float>(particula.gethvx());
  auto hvy    = static_cast<float>(particula.gethvy());
  auto hvz    = static_cast<float>(particula.gethvz());
  auto vx_dat = static_cast<float>(particula.getvx());
  auto vy_dat = static_cast<float>(particula.getvy());
  auto vz_dat = static_cast<float>(particula.getvz());

  return std::make_tuple(px_dat, py_dat, pz_dat, hvx, hvy, hvz, vx_dat, vy_dat, vz_dat);
}

void escribir_datos_particulas(std::ofstream & outputFile, particula const & particula) {
  auto [px_dat, py_dat, pz_dat, hvx, hvy, hvz, vx_dat, vy_dat, vz_dat] = convertirDatos(particula);

  outputFile.write(as_buffer(px_dat), sizeof(px_dat));
  outputFile.write(as_buffer(py_dat), sizeof(py_dat));
  outputFile.write(as_buffer(pz_dat), sizeof(pz_dat));
  outputFile.write(as_buffer(hvx), sizeof(hvx));
  outputFile.write(as_buffer(hvy), sizeof(hvy));
  outputFile.write(as_buffer(hvz), sizeof(hvz));
  outputFile.write(as_buffer(vx_dat), sizeof(vx_dat));
  outputFile.write(as_buffer(vy_dat), sizeof(vy_dat));
  outputFile.write(as_buffer(vz_dat), sizeof(vz_dat));
}


void grid::almacenar_resultados(std::ofstream & outputFile, std::ifstream const & inputFile) {
  // Escribir los parámetros generales
  escribir_parametros_generales(outputFile, inputFile);

  // Escribir los datos de las partículas

  for (int i = 0; i < getnx() * getny() * getnz(); i++) {
    for (auto & particula : bloques[i].particulas) {
        escribir_datos_particulas(outputFile, particula);
    }
  }
}
