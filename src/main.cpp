#include "frame_utils.hpp"
#include <iostream>

int main(){
    Eigen::Vector3d euler_angle;
    euler_angle << 95.0,105.5,260.3;
    euler_angle = euler_angle / 180 * M_PI;
    Eigen::Matrix3d R_b_geodetic = Euler_to_CTM(euler_angle);

    Eigen::Vector3d lat_lon_alt;
    double lat, lon, alt;
    lat = 22.0/180.0*M_PI;
    lon = 153.0/180.0*M_PI;
    alt = 323.0;
    lat_lon_alt << lat,lon,alt;

    Eigen::Vector3d p_eb_e;
    Eigen::Matrix3d R_b_ecef;

    Eigen::Vector3d lat_lon_alt_reconstructed;
    Eigen::Matrix3d R_b_geodetic_reconstructed;

    geodetic2ecef(lat_lon_alt,R_b_geodetic,p_eb_e,R_b_ecef);
    ecef2geodetic(p_eb_e,R_b_ecef,lat_lon_alt_reconstructed,R_b_geodetic_reconstructed);

    std::cout << "lat_lon_alt : \n" << lat_lon_alt << std::endl;
    std::cout << "R_b_geodetic : \n" << R_b_geodetic << std::endl;

    std::cout << "p_eb_e : \n" << p_eb_e << std::endl;
    std::cout << "R_b_ecef : \n" << R_b_ecef << std::endl;

    std::cout << "lat_lon_alt_reconstructed : \n" << lat_lon_alt_reconstructed << std::endl;
    std::cout << "R_b_geodetic_reconstructed : \n" << R_b_geodetic_reconstructed << std::endl;
}