#include <iostream>
#include <fstream>
#include <string>
#include<cstdlib>
#include <vector>
#include <sstream>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <math.h>
#include <ctime>
#include <boost/program_options.hpp>


using std::vector;
using namespace std;
using namespace Eigen;
using namespace boost::program_options;

int main(int argc, char const *argv[])
{   
    boost::program_options::options_description des_cmd("\n Usage: Linear regression based on summary data. \n\n Options: \n");
	des_cmd.add_options()
        ("help,h", "help")
		("file", boost::program_options::value<std::string>(), "Names of the covariance matrices.")
		("phe", boost::program_options::value<std::string>(), "Name of the phenotype that match the column names in the covariace matrix.")
		("covar", boost::program_options::value<std::string>(), "Names of the covariates that match the column names in the covariace matrix, seperated by comma and nothing else. e.g. \"X1,X2,X3\".")
		("out", boost::program_options::value<std::string>(), "prefix of the output file.");

	boost::program_options::variables_map virtual_map;
	try
	{
		boost::program_options::store(boost::program_options::parse_command_line(argc, argv, des_cmd), virtual_map);
	}
	catch (...){ return 0; }
	
	// 定义消息
	boost::program_options::notify(virtual_map);
	
	// 无参数直接返回
	if (virtual_map.empty())
	{
        std::cout << "No parameter is given, exiting..." << std::endl;
		return 0;
	}
	else if (virtual_map.count("help") || virtual_map.count("h"))
	{
		std::cout << des_cmd << std::endl;
		return 0;
	}
	else if (virtual_map.count("file") && virtual_map.count("phe") && virtual_map.count("covar") && virtual_map.count("out"))
	{
        std::cout << "Using parameter: " << std::endl;
        std::cout << "file = " << virtual_map["file"].as<std::string>() << std::endl;
		std::cout << "phe = " << virtual_map["phe"].as<std::string>() << std::endl;
        vector<string> covar = string_split_by_comma(virtual_map["covar"].as<std::string>());
        std::cout << "covar = ";
        for (int i = 0; i < covar.size(); i++)
        {
            std::cout << covar[i];
        }
        
        std::cout << std::endl;
        std::cout << "out = " << virtual_map["out"].as<std::string>() << std::endl;
        std::cout << "Start computing..." << std::endl;
        return 0;
	}
	else
	{
		std::cout << "option error" << std::endl;
        return 0;
	}

    return 0;

    // clock_t start, end;
    // start = clock();
    // string cov_xy_filename = "./test_data/cov_xy_sam.table";
    // string var_x_filename = "./test_data/var_x_sam.table";
    // string cov_yy_filename = "./test_data/cov_yy.table";
    // string result_filename = "./test_data/res.table";
    // vector<string> covar = {"X31.0.0","X1160.0.0", "X1200.0.0", "X1289.0.0",
    //                         "PC1", "PC2", "PC3", "PC4", "PC5"};
    // Matrix<double, Dynamic, Dynamic> Theta;
    // int d = covar.size();
    // Theta.setZero(d+1, d+1);
    // Matrix<double, Dynamic, Dynamic> Theta_0 = Read_matrix_table(cov_yy_filename, covar);

    // Theta(0,0) = Theta_0(0,0);
    // Theta.block(0,2,1,d-1) = Theta_0.block(0,1,1,d-1);
    // Theta.block(2,0,d-1,1) = Theta_0.block(1,0,d-1,1);
    // Theta.block(2,2,d-1,d-1) = Theta_0.block(1,1,d-1,d-1);

    // cal_summary_and_save(cov_xy_filename, var_x_filename, result_filename, Theta, covar);
    // end = clock();
    // double endtime=(double)(end-start)/CLOCKS_PER_SEC;
    // std::cout << "10000 step completed in " << endtime << " seconds." << endl;
    // system("pause");
    // return 0;
}