#include "frame_utils.hpp"


// Function to get the sign of a value
double sign(double value) {
    return (value < 0) ? -1.0 : 1.0;
}


void geodetic2ecef(const Eigen::Vector3d lat_lon_alt, const Eigen::Matrix3d& C_b_geodetic,
                   Eigen::Vector3d& r_eb_e, Eigen::Matrix3d& C_b_ecef){
    
    double lat, lon, alt;
    lat = lat_lon_alt(0);
    lon = lat_lon_alt(1);
    alt = lat_lon_alt(2);

    // Parameters
    double R_0 = 6378137;  // WGS84 Equatorial radius in meters
    double e = 0.0818191908425;  // WGS84 eccentricity

    // Calculate transverse radius of curvature using (2.105)
    double R_E = R_0 / std::sqrt(1 - std::pow(e * std::sin(lat), 2));

    // Convert position using (2.112)
    double cos_lat = std::cos(lat);
    double sin_lat = std::sin(lat);
    double cos_long = std::cos(lon);
    double sin_long = std::sin(lon);
    r_eb_e << (R_E + alt) * cos_lat * cos_long,
              (R_E + alt) * cos_lat * sin_long,
              ((1 - std::pow(e, 2)) * R_E + alt) * sin_lat;

    // Calculate ECEF to NED coordinate transformation matrix using (2.150)
    Eigen::Matrix3d C_ecef_geodetic;
    C_ecef_geodetic << -sin_lat * cos_long, -sin_lat * sin_long, cos_lat,
             -sin_long, cos_long, 0,
             -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;


    // Transform attitude using (2.15)
    C_b_ecef = C_ecef_geodetic.transpose() * C_b_geodetic;
}

void ecef2geodetic(const Eigen::Vector3d& r_eb_e, const Eigen::Matrix3d& C_b_ecef,
                   Eigen::Vector3d& lat_lon_alt, Eigen::Matrix3d& C_b_geodetic){
    // Parameters
    double R_0 = 6378137;  // WGS84 Equatorial radius in meters
    double e = 0.0818191908425;  // WGS84 eccentricity

    // Convert position using Borkowski closed-form exact solution
    // From (2.113) in Paul Groves
    double lon = std::atan2(r_eb_e(1), r_eb_e(0));

    // From (C.29) and (C.30) in Paul Groves
    double k1 = std::sqrt(1 - e * e) * std::abs(r_eb_e(2));
    double k2 = e * e * R_0;
    double beta = std::sqrt(r_eb_e(0) * r_eb_e(0) + r_eb_e(1) * r_eb_e(1));
    double E = (k1 - k2) / beta;
    double F = (k1 + k2) / beta;

    // From (C.31) in Paul Groves
    double P = 4.0 / 3.0 * (E * F + 1);

    // From (C.32) in Paul Groves
    double Q = 2 * (E * E - F * F);

    // From (C.33) in Paul Groves
    double D = std::pow(P, 3) + std::pow(Q, 2);

    // From (C.34) in Paul Groves
    double V = std::pow(std::sqrt(D) - Q, 1.0 / 3.0) - std::pow(std::sqrt(D) + Q, 1.0 / 3.0);

    // From (C.35) in Paul Groves
    double G = 0.5 * (std::sqrt(E * E + V) + E);

    // From (C.36) in Paul Groves
    double T = std::sqrt(G * G + (F - V * G) / (2 * G - E)) - G;

    // From (C.37) in Paul Groves
    double lat = sign(r_eb_e(2)) * std::atan(  (1 - T*T)  /  (2 * T * std::sqrt(1 - e*e))  );

    // From (C.38) in Paul Groves
    double alt = (beta - R_0 * T) * std::cos(lat) +
      ( r_eb_e(2) - sign(r_eb_e(2)) * R_0 * std::sqrt(1 - e*e) ) * std::sin(lat);

    // Convert position using (2.112)
    double cos_lat = std::cos(lat);
    double sin_lat = std::sin(lat);
    double cos_long = std::cos(lon);
    double sin_long = std::sin(lon);

    // Calculate ECEF to NED coordinate transformation matrix using (2.150)
    Eigen::Matrix3d C_ecef_geodetic;
    C_ecef_geodetic << -sin_lat * cos_long, -sin_lat * sin_long, cos_lat,
             -sin_long, cos_long, 0,
             -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;


    // Transform attitude using (2.15)
    C_b_geodetic = C_ecef_geodetic * C_b_ecef;

    lat_lon_alt << lat,lon,alt;
}

Eigen::Matrix3d Euler_to_CTM(const Eigen::Vector3d& eul) {
    double sin_phi = std::sin(eul(0));
    double cos_phi = std::cos(eul(0));
    double sin_theta = std::sin(eul(1));
    double cos_theta = std::cos(eul(1));
    double sin_psi = std::sin(eul(2));
    double cos_psi = std::cos(eul(2));

    Eigen::Matrix3d C;
    C(0, 0) = cos_theta * cos_psi;
    C(0, 1) = cos_theta * sin_psi;
    C(0, 2) = -sin_theta;
    C(1, 0) = -cos_phi * sin_psi + sin_phi * sin_theta * cos_psi;
    C(1, 1) = cos_phi * cos_psi + sin_phi * sin_theta * sin_psi;
    C(1, 2) = sin_phi * cos_theta;
    C(2, 0) = sin_phi * sin_psi + cos_phi * sin_theta * cos_psi;
    C(2, 1) = -sin_phi * cos_psi + cos_phi * sin_theta * sin_psi;
    C(2, 2) = cos_phi * cos_theta;

    return C;
}