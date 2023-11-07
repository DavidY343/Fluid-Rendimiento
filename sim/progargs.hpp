//
// Created by david on 11/7/23.
//

#ifndef ARCOS_PROGARGS_HPP
#define ARCOS_PROGARGS_HPP

#include <array>
#include <fstream>
#include <string>
#include <vector>

namespace constantes {
  double const r_const                   = 1.695;
  int const p_const                      = 1000;
  double const ps_const                  = 3.0;
  int const sc_const                     = 30000;
  double const dv_const                  = 128.0;
  double const u_const                   = 0.4;
  double const dp_const                  = 0.0002;
  double const t_const                   = 0.001;
  std::array<double, 3> const g_const    = {0.0, -9.8, 0.0};
  std::array<double, 3> const bmax_const = {0.065, 0.1, 0.065};
  std::array<double, 3> const bmin_const = {-0.065, -0.08, -0.065};
}  // namespace constantes

class particula {
  public:
    // Constructor
    particula(int id, double px, double py, double pz, double hvx, double hvy, double hvz,
              double vy, double vx, double vz)
      : id(id), px(px), py(py), pz(pz), hvx(hvx), hvy(hvy), hvz(hvz), vy(vy), vx(vx), vz(vz) { }

    // Constructor de copia
    particula(particula const & other) = default;

    // Constructor de movimiento
    particula(particula && other) = default;

    // Destructor
    ~particula() = default;

    // Funciones miembro
    void setPosicion(double new_x, double new_y, double new_z) {
      px = new_x;
      py = new_y;
      pz = new_z;
    }

    [[nodiscard]] double getVelocidadX() const { return vx; }

    [[nodiscard]] double getVelocidadY() const { return vy; }

    [[nodiscard]] double getVelocidadZ() const { return vz; }

    // Operador de asignación de copia
    particula & operator=(particula const & other) {
      if (this != &other) {
        px  = other.px;
        py  = other.py;
        pz  = other.pz;
        hvx = other.hvx;
        hvy = other.hvy;
        hvz = other.hvz;
        vy  = other.vy;
        vx  = other.vx;
        vz  = other.vz;
      }
      return *this;
    }

    // Operador de asignación de movimiento  NO SE
    particula & operator=(particula && other) noexcept {
      if (this != &other) { }
      return *this;
    }

    /*Getters*/
    [[nodiscard]] double getpx() const { return px; }

    [[nodiscard]] double getpy() const { return py; }

    [[nodiscard]] double getpz() const { return pz; }

  private:
    int id;
    double px;
    double py;
    double pz;
    double hvx;
    double hvy;
    double hvz;
    double vy;
    double vx;
    double vz;
};

// funciones
int comprobar_archivos(std::vector<std::string> const & argumentos, std::ifstream const & inputFile,
                       std::ofstream const & outputFile);

int comprobar_params(int argc, std::vector<std::string> const & argumentos,
                     std::ifstream const & inputFile, std::ofstream const & outputFile);

void init_params(std::ifstream & inputFile);

template <typename T>
  requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
char * as_writable_buffer(T & value);

template <typename T>
  requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
char const * as_buffer(T const & value);

template <typename T>
  requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
T read_binary_value(std::istream & is);

#endif  // ARCOS_PROGARGS_H
