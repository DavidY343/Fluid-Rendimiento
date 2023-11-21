//
// Created by david on 11/7/23.
//

#ifndef ARCOS_PROGARGS_HPP
#define ARCOS_PROGARGS_HPP

#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
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
              double vx, double vy, double vz)
      : id(id), px(px), py(py), pz(pz), hvx(hvx), hvy(hvy), hvz(hvz), vx(vx), vy(vy), vz(vz) { }

    // Constructor de copia
    particula(particula const & other) = default;

    // Constructor de movimiento
    particula(particula && other) = default;

    // Destructor
    ~particula() = default;

    // Operador de asignación de copia
    particula & operator=(particula const & other) {
      if (this != &other) {
        id  = other.id;
        px  = other.px;
        py  = other.py;
        pz  = other.pz;
        hvx = other.hvx;
        hvy = other.hvy;
        hvz = other.hvz;
        vx  = other.vx;
        vy  = other.vy;
        vz  = other.vz;
        ax  = other.ax;
        ay  = other.ay;
        az  = other.az;
      }
      return *this;
    }

    // Operador de asignación de movimiento
    particula & operator=(particula && other) noexcept {
      if (this != &other) { }
      return *this;
    }

    void inicializar_densidad_aceleracion();
    void interactuar_densidad(particula & part, double h);

    void interactuar_aceleracion(particula & part, double h, double m);
    void interactuar_aceleracion2(particula & part, std::vector<double> & d_a);
    void transformar_densidad(double h, double m);
    [[nodiscard]] double calcular_factor_aceleracion(double h, double dist, double m,
                                                     particula const & part) const;
    static double calcular_factor_v(double h, double m);
    [[nodiscard]] std::vector<double> calcular_d_a(particula const & part, double h, double dist,
                                                   double m) const;
    void colisionLimiteEjeX(bool lim_inf);
    void colisionLimiteEjeY(bool lim_inf);
    void colisionLimiteEjeZ(bool lim_inf);
    void limiteRecintox(bool lim_inf);
    void limiteRecintoy(bool lim_inf);
    void limiteRecintoz(bool lim_inf);
    [[nodiscard]] double calcularDistancia(particula const & part) const;

    void actualizarMovimiento();

    /*Getters*/
    [[nodiscard]] double getpx() const { return px; }

    [[nodiscard]] double getpy() const { return py; }

    [[nodiscard]] double getpz() const { return pz; }

    [[nodiscard]] double getp() const { return p; }

    [[nodiscard]] int getid() const { return id; }

    [[nodiscard]] double gethvx() const { return hvx; }

    [[nodiscard]] double gethvy() const { return hvy; }

    [[nodiscard]] double gethvz() const { return hvz; }

    [[nodiscard]] double getvx() const { return vx; }

    [[nodiscard]] double getvy() const { return vy; }

    [[nodiscard]] double getvz() const { return vz; }

    [[nodiscard]] double getax() const { return ax; }

    [[nodiscard]] double getay() const { return ay; }

    [[nodiscard]] double getaz() const { return az; }

    void setpx(double new_px) { px = new_px; }

    void setpy(double new_py) { py = new_py; }

    void setpz(double new_pz) { pz = new_pz; }

    void setvx(double new_vx) { vx = new_vx; }

    void setvy(double new_vy) { vy = new_vy; }

    void setvz(double new_vz) { vz = new_vz; }

    void sethvx(double new_hvx) { hvx = new_hvx; }

    void sethvy(double new_hvy) { hvy = new_hvy; }

    void sethvz(double new_hvz) { hvz = new_hvz; }

    void setax(double new_ax) { ax = new_ax; }

    void setay(double new_ay) { ay = new_ay; }

    void setaz(double new_az) { az = new_az; }

    void setp(double new_p) { p = new_p; }

  private:
    int id;
    double px;
    double py;
    double pz;
    double hvx;
    double hvy;
    double hvz;
    double vx;
    double vy;
    double vz;
    double ax{0};
    double ay{0};
    double az{0};
    double p{0};
};

int n_params(int argc);

int comprobar_archivos(std::vector<std::string> const & argumentos, std::ifstream const & inputFile,
                       std::ofstream const & outputFile);

int comprobar_params(std::vector<std::string> const & argumentos, std::ifstream const & inputFile,
                     std::ofstream const & outputFile);

void comprobar_fallos_cabecera(std::vector<particula> const & particulas, int n_particulas_int);
std::vector<particula> crear_particulas(std::ifstream const & inputFile);

template <typename T>
  requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
char * as_writable_buffer(T & value);

template <typename T>
  requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
char const * as_buffer(T const & value);

template <typename T>
  requires(std::is_integral_v<T> or std::is_floating_point_v<T>)
T read_binary_value(std::istream & is);

#endif
