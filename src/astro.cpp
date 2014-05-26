/*Class of astronomy functions*/

#include "astro.hpp"

int Astro::find_bin (double value, double min_value, double bin_size) { 
  //! Find bin corresponding to the input value given the minimum
  //! value in range and bin size.
  return floorf((value - min_value) / bin_size);
}

int Astro::num_bins (double min_value, double max_value, double bin_size) {
  //! Find number of bins in given range for a given bin size.
  return floorf((max_value - min_value) / bin_size);
}

bool Astro::within (double value, double min_value, double max_value) {
  //! Deterimine whether or not a values is within the limits 
  //! provided.
  return (value >= min_value && value < max_value);
}

double Astro::deg2rad (double angle) {
  //! Function that converts angle from degrees to radians.
  return angle * M_PI / 180.0;
}

double Astro::rad2deg (double angle) {
  //! Function that converts angle from radians to degrees.
  return angle * 180.0 / M_PI;
}

double Astro::angsep (double ra1, double dec1, double ra2, double dec2) {
  //! Function that returns the angular separation (in radians) between two points.
  if(ra1 == ra2 && dec1 == dec2)
    return 0.0;
  else
    return acos(sin(deg2rad(dec1)) * sin(deg2rad(dec2)) + cos(deg2rad(dec1)) * 
		cos(deg2rad(dec2)) * cos(deg2rad(ra1) - deg2rad(ra2)));
}

double Astro::mean (const std::vector<double> &elements) {
  //! Function that computes the mean value of a vector of doubles.
  double sum = std::accumulate(elements.begin(), elements.end(), 0.0);
  return sum / double(elements.size());
}

double Astro::median (std::vector<double> elements) {
  //! Function that computes the median value of a vector of doubles.
  //! Pass vector by value to leave original vector unaltered.
  int size = elements.size();
  double median;
  std::sort(elements.begin(), elements.end());
  if (size % 2 == 0)
    median = (elements[size / 2 - 1] + elements[size / 2]) / 2;
  else
    median = elements[size / 2];
  return median;
}

double Astro::min (const std::vector<double> &elements) {
  //! Function that computes the minimum value of a vector of doubles.
  return *std::min_element(elements.begin(), elements.end());
}

double Astro::max (const std::vector<double> &elements) {
  //! Function that computes the maximum value of a vector of doubles.
  return *std::max_element(elements.begin(), elements.end());
}