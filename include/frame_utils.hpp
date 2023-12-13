#include <cmath>
#include <Eigen/Dense> 
#include <iostream>

void geodetic2ecef(const Eigen::Vector3d lat_lon_alt, const Eigen::Matrix3d& C_b_geodetic,
                   Eigen::Vector3d& r_eb_e, Eigen::Matrix3d& C_b_ecef);

void ecef2geodetic(const Eigen::Vector3d& r_eb_e, const Eigen::Matrix3d& C_b_ecef,
                   Eigen::Vector3d& lat_lon_alt, Eigen::Matrix3d& C_b_geodetic);

Eigen::Matrix3d Euler_to_CTM(const Eigen::Vector3d& eul);
