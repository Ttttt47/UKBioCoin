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
#include "tools/cpp_cdfs-master/cdf_chisqt/cdf_base.h"

#include <boost/program_options.hpp>


using std::vector;
using namespace std;
using namespace Eigen;
using namespace boost::program_options;


Vector<double, Dynamic> cal_summary(Matrix<double, Dynamic, Dynamic> Theta)
{   
    int n = 270000;  // estimated. effective sample size.
    int d = Theta.rows();
    Matrix<double, Dynamic, Dynamic> Omega = Theta.block(1,1,d-1,d-1);
    Matrix<double, Dynamic, Dynamic> Lambda = Omega.diagonal().asDiagonal();
    Matrix<double, Dynamic, Dynamic> b = Theta.block(1,0,d-1,1).array() / Omega.diagonal().array();
    
    Matrix<double, Dynamic, Dynamic> rv_Omega = Omega.inverse();
    Matrix<double, Dynamic, Dynamic> beta = rv_Omega * Lambda * b;
    Matrix<double, Dynamic, Dynamic> temp = (Theta.block(0,0,1,1) - (beta.transpose() * Lambda * b))/(n-(d-2)-1);
    Matrix<double, Dynamic, Dynamic> var_beta = temp(0,0) * rv_Omega;

    // beta and SE.
    Vector<double, Dynamic> t_stats = beta.array() / var_beta.diagonal().array().sqrt().array();
    Vector<double, Dynamic> log_p_val;
    Vector<double, Dynamic> p_val;
    Vector<double, Dynamic> log10_p_val;

    log_p_val.setZero(t_stats.size());
    log10_p_val = log_p_val;
    for (int i = 0; i < t_stats.size(); i++)
    {   
        log_p_val[i] = cdf_t_log(-abs(t_stats[i]), n-2) + log(2);
        log10_p_val[i] = log_p_val[i] * log10(exp(1));
    }
    

    // std::cout << rv_Omega << endl;
    // std::cout << t_stats << endl;
    // std::cout << log10_p_val << endl;
    // .asDiagonal();
    Vector<double, 4> summary = {beta(0,0), sqrt(var_beta.diagonal()(0,0)),
                                 t_stats[0], log10_p_val[0]};
    return(summary);
}

Matrix<double, Dynamic, Dynamic> Read_matrix_table(string filename, vector<string> colnames)
{
    ifstream fin;
    fin.open(filename, ios::in);
    if( !fin )
    {   
        cout << "Error opening " << filename << " for input" << endl;   
        exit(-1);  
    }
    string line;
    string word;
    double el=0;
    vector<double> line_data;
    int length_covar=colnames.size();
    Vector<int, Dynamic> covar_id; 
    covar_id = covar_id.setOnes(length_covar)*-1;

    Matrix<double, Dynamic, Dynamic> Theta;
    
    getline(fin, line);  // skip first line.
    istringstream iline(line);
    int word_id = 0;
    while (iline >> word)
    {   
        word = word.substr(1,word.size()-2);
        bool flag = false;
        for (int i = 0; i < length_covar; i++)
        {
            if (word == colnames[i])
            {   
                covar_id[i] = word_id;
                break;
            }
        }
        word_id++;
    }
    // std::cout << covar_id << endl;
    
    Theta.setZero(length_covar, length_covar);
    int line_id = 0;
    while(getline(fin, line))
    {   
        istringstream iline(line);
        iline >> word;  // ignore first colunm.
        word_id = 0;
        while (iline >> word)
        {   
            for (int j = 0; j < covar_id.size(); j++)
            {
                if (line_id == covar_id[j])
                {
                    for (int i = 0; i < covar_id.size(); i++)
                    {
                        if (word_id == covar_id[i])
                        {
                            if (word == "NA")
                            {
                                el = std::numeric_limits<double>::quiet_NaN();
                            }
                            else
                            {
                                el = std::stod(word);
                            }
                            Theta(j,i) = el;
                            break;
                        }
                    }
                }
            }
            

            word_id ++;
        }
        line_id ++;
    }
    // std::cout << Theta.col(1) << endl;
    // std::cout << Theta.row(1) << endl;
    return(Theta);
    // exit(0);
}

void cal_summary_and_save(string cov_xy_filename, string var_x_filename, string result_filename, Matrix<double, Dynamic, Dynamic> Theta, vector<string> covar)
{
    ofstream fout;
    fout.open(result_filename, ios::out);
    if( !fout )
    {   
        cout << "Error opening " << result_filename << " for output" << endl;   
        exit(-1);  
    }
    ifstream fin;
    fin.open(cov_xy_filename, ios::in);
    if( !fin )
    {   
        cout << "Error opening " << cov_xy_filename << " for input" << endl;   
        exit(-1);  
    }
    ifstream fin_x;
    fin_x.open(var_x_filename, ios::in);
    if( !fin_x )
    {   
        cout << "Error opening " << var_x_filename << " for input" << endl;   
        exit(-1);  
    }
    string line;
    string line_x;
    string word;
    string word_x;
    double el=0;
    vector<double> line_data;
    
    int length_covar=covar.size();
    Vector<int, Dynamic> covar_id; 
    covar_id = covar_id.setOnes(length_covar)*-1;
    getline(fin, line); 
    getline(fin_x, line_x);

    istringstream iline(line);
    int word_id = 0;
    while (iline >> word)
    {   
        word = word.substr(1,word.size()-2);
        bool flag = false;
        for (int i = 0; i < length_covar; i++)
        {
            if (word == covar[i])
            {   
                covar_id[i] = word_id;
                break;
            }
        }
        word_id++;
    }

    int line_id = 0;
    fout << "BETA" << ' ' << "SE" << ' ' << "T-STAT" << ' ' << "log10_P" << '\n';
    while(getline(fin, line) && getline(fin_x, line_x))
    {   
        istringstream iline(line);
        iline >> word;  // ignore first colunm.
        istringstream iline_x(line_x);
        iline_x >> word_x;  // ignore first colunm.

        word_id = 0;
        // std::cout << Theta.row(1) << endl;
        while (iline >> word)
        {   
            for (int i = 0; i < covar_id.size(); i++)
            {
                if (word_id == covar_id[i])
                {
                    if (word == "NA")
                    {
                        el = std::numeric_limits<double>::quiet_NaN();
                    }
                    else
                    {
                        el = std::stod(word);
                    }
                    Theta(1, i>0 ? i+1 : 0) = el;
                    break;
                }
            }
            word_id ++;
        }
        Theta.col(1) = Theta.row(1);
        iline_x >> word_x;
        if (word_x == "NA")
        {
            el = std::numeric_limits<double>::quiet_NaN();
        }
        else
        {
            el = std::stod(word_x);
        }
        Theta(1,1) = el;
        // std::cout << Theta.row(0) << endl;
        // std::cout << Theta.row(1) << endl;
        // calculating summary data 
        Vector<double, Dynamic> s;
        Vector<double, 4> summary = cal_summary(Theta);
        fout << summary[0] << ' ' << summary[1] << ' ' << summary[2] << ' ' << summary[3] << '\n';
        // writing results to table.
        // std::cout << summary << endl;
        line_id ++;
        if (line_id%1000==0)
        {
            std::cout << line_id << endl;
        }
    }
    fout.close();
}

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
            std::cout << covar[i] << " ";
        }
        
        std::cout << std::endl;
        std::cout << "out = " << virtual_map["out"].as<std::string>() << std::endl;
        std::cout << "Start computing..." << std::endl;
	}
	else
	{
		std::cout << "option error" << std::endl;
        return 0;
	}

    clock_t start, end;
    start = clock();
    string cov_xy_filename = virtual_map["file"].as<std::string>() + "_cov_xy.table";
    string var_x_filename = virtual_map["file"].as<std::string>() + "_var_x.table";
    string cov_yy_filename = virtual_map["file"].as<std::string>() + "_cov_yy.table";
    string result_filename = virtual_map["out"].as<std::string>() + "_results.table";
    vector<string> covar = string_split_by_comma(virtual_map["covar"].as<std::string>());
    // vector<string> covar = {"X31.0.0","X1160.0.0", "X1200.0.0", "X1289.0.0",
    //                         "PC1", "PC2", "PC3", "PC4", "PC5"};
    Matrix<double, Dynamic, Dynamic> Theta;
    int d = covar.size();
    Theta.setZero(d+1, d+1);
    Matrix<double, Dynamic, Dynamic> Theta_0 = Read_matrix_table(cov_yy_filename, covar);

    Theta(0,0) = Theta_0(0,0);
    Theta.block(0,2,1,d-1) = Theta_0.block(0,1,1,d-1);
    Theta.block(2,0,d-1,1) = Theta_0.block(1,0,d-1,1);
    Theta.block(2,2,d-1,d-1) = Theta_0.block(1,1,d-1,d-1);

    cal_summary_and_save(cov_xy_filename, var_x_filename, result_filename, Theta, covar);
    end = clock();
    double endtime=(double)(end-start)/CLOCKS_PER_SEC;
    std::cout << "Completed in " << endtime << " seconds." << endl;

    return 0;
}

// .\main.exe --file test_data/sam --phe X31.0.0 --covar X1160.0.0,X1200.0.0,X1289.0.0,PC1,PC2,PC3,PC4,PC5 --out test_data/test