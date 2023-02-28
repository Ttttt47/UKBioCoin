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

vector<string> string_split_by_comma(string s){
	int n = s.size();
	for (int i = 0; i < n; ++i){
		if (s[i] == ','){
			s[i] = ' ';
		}
	}
	istringstream out(s);
    string str;
	vector<string> strs;
	while (out >> str){
		strs.push_back(str);
	}
    return(strs);
}


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
}