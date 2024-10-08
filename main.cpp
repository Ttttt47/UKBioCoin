#include <iostream>
#include <fstream>
#include <string>
#include<cstdlib>
#include <vector>
#include <sstream>
#include <limits>
#include <Eigen/Dense>
#include "tools/cpp_cdfs-master/cdf_chisqt/cdf_base.h"
#include <Eigen/Core>
#include <math.h>
#include <ctime>
#include <boost/program_options.hpp>

using std::vector;
using namespace std;
using namespace Eigen;
using namespace boost::program_options;



double cal_VIF(Matrix<double, Dynamic, Dynamic> Omega)
{   
    // cal_VIF is used to calculate the variance inflation factor of the estimated parameters.
    // Omega is the estimated covariance matrix of the parameters.
    int d = Omega.rows();
    double R_2 = (Omega.row(0).segment(1,d-1) * Omega.block(1,1,d-1,d-1).inverse() * Omega.col(0).segment(1,d-1)).value()/Omega(0,0);
    double VIF = 1/(1-R_2);
    return(VIF);
}

Vector<double, Dynamic> cal_summary(Matrix<double, Dynamic, Dynamic> Theta, Matrix<double, Dynamic, Dynamic> D_inv, int n)
{   
    int d = Theta.rows();
    Matrix<double, Dynamic, Dynamic> Omega = Theta.block(1,1,d-1,d-1);
    Matrix<double, Dynamic, Dynamic> Omega_inv = Omega;
    // calculate the inverse of Omega by inverting the truncked matrix & D_inv.
    double q = 1/(Omega(0,0) - (Omega.block(0,1,1,d-2) * D_inv * Omega.block(1,0,d-2,1)).value());
    Omega_inv(0,0) = q;
    Omega_inv.block(0,1,1,d-2) = -q * Omega.block(0,1,1,d-2) * D_inv;
    Omega_inv.block(1,0,d-2,1) = -q * D_inv * Omega.block(1,0,d-2,1);
    Omega_inv.block(1,1,d-2,d-2) =  D_inv + q * D_inv * Omega.block(1,0,d-2,1) * Omega.block(0,1,1,d-2) * D_inv;

    Vector<double, Dynamic> a = Theta.block(0,1,1,d-1).transpose();
    Vector<double, Dynamic> theta = Omega_inv * a;

    double rss = n * (1 - (theta.transpose() * a).value())/(n-(d-2)-1);
    double var_beta = rss * Omega_inv(0,0) / n;
    double t_stat = theta[0] / sqrt(var_beta);
    double log_p_val = cdf_t_log(-abs(t_stat), n-2) + log(2);
    double log10_p_val = log_p_val * log10(exp(1));
    double VIF = var_beta / (rss/((n-1)*Theta(1,1)));

    Vector<double, 5> summary = {theta[0], sqrt(var_beta), t_stat, -log10_p_val, VIF};
    return(summary);
}

Matrix<double, Dynamic, Dynamic> Read_matrix_table(string filename, vector<string> colnames, string y_missing_file,bool use_missing_rate_estim, double &non_missing_rate)
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
    
    getline(fin, line);  // colname lines
    istringstream iline(line);
    int word_id = 0;
    while (iline >> word)
    {   
        word = word.substr(1,word.size()-2);  // remove " and " around the covar name.
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
    bool err_flag = false;
    string cov_notfind = "";
    for (int i = 0; i < length_covar; i++){
        if (covar_id[i] == -1){
            err_flag = true;
            cov_notfind = cov_notfind + " " + colnames[i];
        }
    }
    if (err_flag){
        cout << "Error: phenotypes/covariates" << cov_notfind << " not find in the input file." << endl;   
        exit(-1); 
    }

    if (y_missing_file!="" && use_missing_rate_estim){
        ifstream fin_miss;
        fin_miss.open(y_missing_file, ios::in);
        if( !fin_miss )
        {   
            cout << "Error opening " << y_missing_file << " for input" << endl;   
            exit(-1);  
        }
        string line_miss;
        string word;

        Matrix<double, Dynamic, Dynamic> Theta;
        
        while(getline(fin_miss, line_miss))
        {

            istringstream iline(line_miss);
            iline >> word;
            word = word.substr(1,word.size()-2);
            for (int i = 0; i < length_covar; i++)
            {   
                if (word == colnames[i])
                {   
                    iline >> word;
                    non_missing_rate = non_missing_rate*(1-std::stod(word));
                    break;
                }
            }
        }
    }

    
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
    return(Theta);
}

void cal_summary_and_save(string cov_xy_filename, string var_x_filename, string meta_filename, string result_filename,string x_missing_filename,string y_missing_filename ,Matrix<double, Dynamic, Dynamic> Theta,
                             vector<string> covar, int n, double c, double quality_score, bool use_missing_rate_estim)
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

    ifstream fin_meta;
    fin_meta.open(meta_filename, ios::in);
    if( !fin_meta )
    {   
        cout << "Error opening " << meta_filename << " for input" << endl;   
        exit(-1);  
    }
    ifstream fin_missing_x;
    
    if (use_missing_rate_estim){
        fin_missing_x.open(x_missing_filename, ios::in);
        if( !fin_missing_x )
        {   
            cout << "Error opening " << x_missing_filename << " for input" << endl;   
            exit(-1);  
        }
    }

    string line;
    string line_x;
    string line_meta;
    string line_miss;
    string word;
    string word_x;
    string word_meta;
    string word_miss;
    double el=0;
    double c1;
    vector<double> line_data;
    
    int length_covar=covar.size();
    Vector<int, Dynamic> covar_id; 
    covar_id = covar_id.setOnes(length_covar)*-1;
    getline(fin, line); 
    getline(fin_x, line_x);
    getline(fin_meta, line_meta);
    getline(fin_missing_x, line_miss);

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
    istringstream iline_meta(line_meta);
    while (iline_meta >> word_meta)
    {
        fout << word_meta << ' ';
    }
    fout << "BETA" << ' ' << "SE" << ' ' << "T-STAT" << ' ' << "-log10_P" << ' ' << "VIF" << ' ' << "nobs" << ' ' << "Quality-Score" << '\n';
    while(getline(fin, line) && getline(fin_x, line_x) && getline(fin_meta, line_meta))
    {   
        
        if (use_missing_rate_estim){
            getline(fin_missing_x, line_miss);
            istringstream iline_miss(line_miss);
            iline_miss >> word_miss;  // ignore first colunm.
            iline_miss >> word_miss; 
            c1 = std::stod(word_miss);
        }
        istringstream iline_meta(line_meta);
        while (iline_meta >> word_meta)
        {
            fout << word_meta << ' ';
        }
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
        int size;
        if (use_missing_rate_estim)
        {   
            size = floor(n*quality_score*(1-c1));

        }else
        {
            size = floor(n*c);
        }
        int d=Theta.rows();

        Matrix<double, Dynamic, Dynamic> D_inv = Theta.block(2,2,d-2,d-2).inverse();
        // check if D is invertible. if not, then the calculation will be halt.
        if (D_inv.hasNaN())
        {
            std::cout << "Sigularity detected in the covariance matrix, please check if the covariates are linearly dependent." << endl;
            exit(-1);
        }
        Vector<double, 5> summary = cal_summary(Theta, D_inv, size);
        
        // double VIF = cal_VIF(Theta.block(1,1,d-1,d-1)); // depriciated methods for calculating VIF.
        // double VIF = 0;
        if (use_missing_rate_estim)
        {
            fout << summary[0] << ' ' << summary[1] << ' ' << summary[2] << ' ' << summary[3] << ' ' << summary[4] << ' ' << ' ' << size << ' ' << quality_score*(1-c1) << '\n';
        }else
        {
            fout << summary[0] << ' ' << summary[1] << ' ' << summary[2] << ' ' << summary[3] << ' ' << summary[4] << ' ' << ' ' << size << ' ' << c << '\n';
        }
        
        // writing results to table.
        // std::cout << summary << endl;
        line_id ++;
        if (line_id%10000==0)
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
    std::cout << "-----------------------------------" << std::endl;
    std::cout << "UKBioCoin (UKC) engine Version: V1.3" << std::endl;
    std::cout << "Written by: Jing-cheng He" << std::endl;
    std::cout << "GitHub: https://github.com/Ttttt47/UKBioCoin" << std::endl;
    std::cout << "Email: jc_he@zju.edu.cn" << std::endl;
    std::cout << "-----------------------------------" << std::endl;
    
    boost::program_options::options_description des_cmd("\n Usage: Linear regression based on summary data. \n\n Options: \n");
	des_cmd.add_options()
        ("help,h", "help")
        ("use-missing-rate-file", "Whether missing rate files are available.")
        ("use-missing-rate-estimate", "Whether use missing_rate files to estimate sample size.")
		("file", boost::program_options::value<std::string>(), "Names of the covariance matrices.")
		("phe", boost::program_options::value<std::string>(), "Name of the phenotype that match the column names in the covariace matrix.")
		("covar", boost::program_options::value<std::string>(), "Names of the covariates that match the column names in the covariace matrix, seperated by comma and nothing else. e.g. \"X1,X2,X3\".")
		("out", boost::program_options::value<std::string>(), "prefix of the output file.")
        ("totalsize", boost::program_options::value<int>(), "Total Sample size of the regression, default is 292216.")
        ("overall-non-missing-rate", boost::program_options::value<double>(), "Overall non-missing rate, default is 0.9.");
	boost::program_options::variables_map virtual_map;
	try
	{
		boost::program_options::store(boost::program_options::parse_command_line(argc, argv, des_cmd), virtual_map);
	}
	catch (...){ return 0; }

	boost::program_options::notify(virtual_map);
	
    int n;
    double c=1;
    bool use_missing_rate_estim;

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
	else if (virtual_map.count("file") && virtual_map.count("phe") && virtual_map.count("out"))
	{
        std::cout << "Using parameter: " << std::endl;
        std::cout << "file = " << virtual_map["file"].as<std::string>() << std::endl;
		std::cout << "phe = " << virtual_map["phe"].as<std::string>() << std::endl;
        if (virtual_map.count("covar"))
        {
            vector<string> covar = string_split_by_comma(virtual_map["covar"].as<std::string>());
            std::cout << "covar = ";
            for (int i = 0; i < covar.size(); i++)
            {
                std::cout << covar[i] << " ";
                
            }
            std::cout << std::endl;
        }else
        {
            std::cout << "No covar is given." << std::endl;
        }

        if (virtual_map.count("totalsize"))
        {
            n = virtual_map["totalsize"].as<int>();
        }
        else{
            n = 292216;  // default total sample size for UKB.
        }

        use_missing_rate_estim = virtual_map.count("use-missing-rate-estimate");

        if (use_missing_rate_estim)
        {
            std::cout << "using estimated Sample size." << std::endl;
        } else
        {
            if (virtual_map.count("overall-non-missing-rate"))
            {
                c = virtual_map["overall-non-missing-rate"].as<double>();
            }
            else{
                c = 0.9;  // default overall missing rate.
            }
            std::cout << "Sample size set as " << floor(c*n) << std::endl;
        }
        
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
    string meta_filename = virtual_map["file"].as<std::string>() + "_meta.table";
    string result_filename = virtual_map["out"].as<std::string>() + "_results.table";
    string x_missing_filename = "";
    string y_missing_filename = "";
    vector<string> covar;
    if (virtual_map.count("covar"))
    {
        covar = string_split_by_comma(virtual_map["phe"].as<std::string>()+","+virtual_map["covar"].as<std::string>());
    }else{
        covar = string_split_by_comma(virtual_map["phe"].as<std::string>());
    }

    if (use_missing_rate_estim)
    {
        x_missing_filename = virtual_map["file"].as<std::string>() + "_x_missing.table";
        y_missing_filename = virtual_map["file"].as<std::string>() + "_y_missing.table";
    }

    Matrix<double, Dynamic, Dynamic> Theta;
    int d = covar.size();
    double quality_score = 1;
    Theta.setZero(d+1, d+1);
    Matrix<double, Dynamic, Dynamic> Theta_0 = Read_matrix_table(cov_yy_filename, covar, y_missing_filename, use_missing_rate_estim, quality_score);


    Theta(0,0) = Theta_0(0,0);
    Theta.block(0,2,1,d-1) = Theta_0.block(0,1,1,d-1);
    Theta.block(2,0,d-1,1) = Theta_0.block(1,0,d-1,1);
    Theta.block(2,2,d-1,d-1) = Theta_0.block(1,1,d-1,d-1);

    cal_summary_and_save(cov_xy_filename, var_x_filename, meta_filename, result_filename, x_missing_filename, y_missing_filename, Theta, covar, n, c, quality_score, use_missing_rate_estim);
    end = clock();
    double endtime=(double)(end-start)/CLOCKS_PER_SEC;
    std::cout << "Completed in " << endtime << " seconds." << endl;

    return 0;
}

// .\main.exe --file test_data/sam_10M --phe X31.0.0 --covar X1160.0.0,X1200.0.0,X1289.0.0,PC1,PC2,PC3,PC4,PC5 --out test_data/test --overall-non-missing-rate 1 --use-missing-rate-file --use-missing-rate-estimate
// .\main.exe --file UKB_R/136_1k --phe X102.0.0 --covar X21001.0.0,X1488.0.0,X24021.0.0,PC1,PC2,PC3,PC4,PC5 --out UKB_R/C_imp_res