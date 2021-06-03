/*Class for storing cluster properties*/

#include <algorithm>
#include <vector>
#include "cluster_class.hpp"

void Cluster::add_gal (Galaxy *gal) {
  // Add galaxy properties to Cluster instance.
  mem.push_back(gal);
}

void Cluster::assign_dist (double c, double H0, double Omega_M, double Omega_L) {
  // Calculate angular diameter distance for Cluster instance;
  if (c <= 0)
    throw BadArgumentException("Cluster::assign_dist", "c", "> 0.0");
  if (H0 <= 0)
    throw BadArgumentException("Cluster::assign_dist", "H0", "> 0.0");
  if (z <= 0)
    throw RuntimeException("Cluster::assign_dist", "z", "> 0.0");
  cosmo.set_up(Omega_M, Omega_L);
  da = ((c / H0) * cosmo.angdidis(z));
}

void Cluster::assign_props () {
  // Assign properties to Cluster instance.
  double sum = 0.0;
  std::vector<double> g_ra, g_dec, g_z;
  ngal = mem.size();
  if(ngal > 0) {
    for(int i = 0; i < ngal; i++) {
      g_ra.push_back(mem[i]->P.P[0]);
      g_dec.push_back(mem[i]->P.P[1]);
      g_z.push_back(mem[i]->z);
    }
    ra = astro.median(g_ra);
    dec = astro.median(g_dec);
    z = astro.median(g_z);
    ra_err = astro.stderr_median(g_ra);
    dec_err = astro.stderr_median(g_dec);
    z_err = astro.stderr_median(g_z);
    for(int i = 0; i < ngal; i++)
      sum += astro.rad2deg(astro.angsep(ra, dec, g_ra[i], g_dec[i])); //此处将rad转化为了deg
    size = (sum * 60.0) / double(ngal); //因为默认单位为arcmin，所以乘以60，原本为deg单位
    area = M_PI * pow(size, 2);
  }
}

double Cluster::assign_sn (double bg_expect) {
  // Assign singal-to-noise to Cluster instance
  // given the expected background counts.
  if (bg_expect < 0)
    throw BadArgumentException("Cluster::assign_sn", "bg_expect", ">= 0.0");
  if (ngal <= 0)
    throw RuntimeException("Cluster::assign_sn", "ngal", "> 0.0");
  if (area < 0)
    throw RuntimeException("Cluster::assign_sn", "area", ">= 0.0");
  if(bg_expect == 0)
    sn = -1.0;
  else
    sn = (double(ngal) - (area * bg_expect)) / pow((area * bg_expect), 0.5);

    return sn;
}

void Cluster::update_size (const std::string size_units) {
  // Change cluster size units.
  if (size_units == "deg")
    size = size / 60.0;
  else if (size_units == "Mpc")
    size = astro.deg2rad(size / 60.0) * da;
  area = M_PI * pow(size, 2);
}

void Cluster::clear () {
  // Clear all members from Cluster instance.
  mem.clear();
}

void Cluster::rename (int num_val) {
  // Reset Cluster instance number.
  num = num_val;
}

void Cluster::unique () {
  // Remove duplicate elements from a cluster instance.
  std::sort(mem.begin(), mem.end());
  mem.erase(std::unique(mem.begin(), mem.end()), mem.end());
}
