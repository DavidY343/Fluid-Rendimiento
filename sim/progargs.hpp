//
// Created by david on 11/7/23.
//

#ifndef ARCOS_PROGARGS_HPP
#define ARCOS_PROGARGS_HPP

#include <array>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <iostream>
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
      : id(id), px(px), py(py), pz(pz), hvx(hvx), hvy(hvy), hvz(hvz), vy(vy), vx(vx), vz(vz) {
      ax = 0;
      ay = 0;
      az = 0;
      p = 0;
    }

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
        ax = other.ax;
        ay = other.ay;
        az = other.az;
      }
      return *this;
    }

    // Operador de asignación de movimiento  NO SE
    particula & operator=(particula && other) noexcept {
      if (this != &other) { }
      return *this;
    }

    void inicializar_densidad_aceleracion(){
      p = 0;
      ax = constantes::g_const[0];
      ay = constantes::g_const[1];
      az = constantes::g_const[2];
    }

    void imprimir_datos(){
      using namespace std;
      cout<<"particula con id:"<<id<<" en posicion:("<<px<<", "<<py<<", "<<pz<<")"<<endl;
    }

    void interactuar_densidad(particula part, double h, bool sumar_a_ambas_part){  //hay q hacerlo sin tener q pasar h como constante, q pereza
      double distancia_cuadrado = pow(this->px - part.px, 2) + pow(this->py - part.py, 2) + pow(this->pz - part.pz, 2);
      if(distancia_cuadrado< pow(h, 2)){
        p += pow(pow(h, 2) - distancia_cuadrado, 3);
      }
      if(sumar_a_ambas_part){
        part.p += pow(pow(h, 2) - distancia_cuadrado, 3);
      }
    }

    void interactuar_aceleracion(particula part, double h, double m, bool sumar_a_ambas_part){  //hay q hacerlo sin tener q pasar h como constante, q pereza
      double d = pow(std::max(pow(this->px - part.px, 2) + pow(this->py - part.py, 2) + pow(this->pz - part.pz, 2),pow(10,-12)), 0.5);
      std::vector<double> d_a;
      d_a.push_back((((px-part.px)*((15*m*pow(h-d, 2)*(p + part.p - constantes::p_const))/(M_PI* pow(h, 6)*d)))+(part.vx-vx)*(45/(M_PI*pow(h, 6)*m*constantes::u_const)))/(p*part.p));
      d_a.push_back((((py-part.py)*((15*m*pow(h-d, 2)*(p + part.p - constantes::p_const))/(M_PI* pow(h, 6)*d)))+(part.vy-vy)*(45/(M_PI*pow(h, 6)*m*constantes::u_const)))/(p*part.p));
      d_a.push_back((((pz-part.pz)*((15*m*pow(h-d, 2)*(p + part.p - constantes::p_const))/(M_PI* pow(h, 6)*d)))+(part.vz-vz)*(45/(M_PI*pow(h, 6)*m*constantes::u_const)))/(p*part.p));
      ax += d_a[0];
      ay += d_a[1];
      az += d_a[2];
      if(sumar_a_ambas_part) {
        part.ax -= d_a[0];
        part.ay -= d_a[1];
        part.az -= d_a[2];
      }
    }

    void transformar_densidad(double h){
      p = (p + pow(h, 6)) *(315/(64 * M_PI * pow(h, 9)));
    }

    void colisionLimiteEjeX(bool lim_inf);
    void colisionLimiteEjeY(bool lim_inf);
    void colisionLimiteEjeZ(bool lim_inf);
    void limiteRecintox(bool lim_inf);
    void limiteRecintoy(bool lim_inf);
    void limiteRecintoz(bool lim_inf);

    void actualizarMovimiento();

    /*Getters*/
    [[nodiscard]] double getpx() const { return px; }

    [[nodiscard]] double getpy() const { return py; }

    [[nodiscard]] double getpz() const { return pz; }

    [[nodiscard]] double getp() const { return p; }

    [[nodiscard]] double getid() const { return id; }

    [[nodiscard]] double gethvx() const { return hvx; }

    [[nodiscard]] double gethvy() const { return hvy; }

    [[nodiscard]] double gethvz() const { return hvz; }

    [[nodiscard]] double getvx() const { return vx; }

    [[nodiscard]] double getvy() const { return vy; }

    [[nodiscard]] double getvz() const { return vz; }

    [[nodiscard]] double getax() const { return ax; }

    [[nodiscard]] double getay() const { return ay; }

    [[nodiscard]] double getaz() const { return az; }

    void setpx(double new_ax) { ax = new_ax; }

    void setpy(double new_ay) { ay = new_ay; }

    void setpz(double new_az) { az = new_az; }

    void setvx(double new_ax) { ax = new_ax; }

    void setvy(double new_ay) { ay = new_ay; }

    void setvz(double new_az) { az = new_az; }

    void sethvx(double new_ax) { ax = new_ax; }

    void sethvy(double new_ay) { ay = new_ay; }

    void sethvz(double new_az) { az = new_az; }

    void setax(double new_ax) { ax = new_ax; }

    void setay(double new_ay) { ay = new_ay; }

    void setaz(double new_az) { az = new_az; }


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
    double ax;
    double ay;
    double az;
    double p;
};

// funciones
int comprobar_archivos(std::vector<std::string> const & argumentos, std::ifstream const & inputFile,
                       std::ofstream const & outputFile);

int comprobar_params(int argc, std::vector<std::string> const & argumentos,
                     std::ifstream const & inputFile, std::ofstream const & outputFile);

std::vector<double> longitud_masa(std::ifstream const & inputFile);
int comprobar_fallos_cabecera(std::vector<particula> const & particulas, int n_particulas_int);
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

void escribir_parametros_generales(std::ofstream & outputFile, int num_particulas);


#endif
