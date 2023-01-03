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

using std::vector;
using namespace std;
using namespace Eigen;


//输出空行
void OutPutAnEmptyLine()
{
    cout<<"\n";
}
 
//读取方式: 逐词读取, 词之间用空格区分
//read data from the file, Word By Word
//when used in this manner, we'll get space-delimited bits of text from the file
//but all of the whitespace that separated words (including newlines) was lost. 
void ReadDataFromFileWBW()
{
    ifstream fin("data.txt");  
    string s;  
    while( fin >> s ) 
    {    
        cout << "Read from file: " << s << endl;  
    }
}
 
//读取方式: 逐行读取, 将行读入字符数组, 行之间用回车换行区分
//If we were interested in preserving whitespace, 
//we could read the file in Line-By-Line using the I/O getline() function.
void ReadDataFromFileLBLIntoCharArray()
{
    ifstream fin("data.txt"); 
    const int LINE_LENGTH = 100; 
    char str[LINE_LENGTH];  
    while( fin.getline(str,LINE_LENGTH) )
    {    
        cout << "Read from file: " << str << endl;
    }
}
 
//读取方式: 逐行读取, 将行读入字符串, 行之间用回车换行区分
//If you want to avoid reading into character arrays, 
//you can use the C++ string getline() function to read lines into strings
void ReadDataFromFileLBLIntoString()
{
    ifstream fin("data.txt");  
    string s;  
    while( getline(fin,s) )
    {    
        cout << "Read from file: " << s << endl; 
    }
}
 
//带错误检测的读取方式
//Simply evaluating an I/O object in a boolean context will return false 
//if any errors have occurred
void ReadDataWithErrChecking()
{
    string filename = "dataFUNNY.txt";  
    ifstream fin( filename.c_str());  
    if( !fin ) 
    {   
        cout << "Error opening " << filename << " for input" << endl;   
        exit(-1);  
    }
}

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


int main()
{   
    clock_t start, end;
    start = clock();
    string cov_xy_filename = "./test_data/cov_xy_sam.table";
    string var_x_filename = "./test_data/var_x_sam.table";
    string cov_yy_filename = "./test_data/cov_yy.table";
    string result_filename = "./test_data/res.table";
    vector<string> covar = {"X31.0.0","X1160.0.0", "X1200.0.0", "X1289.0.0",
                            "PC1", "PC2", "PC3", "PC4", "PC5"};
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
    std::cout << "10000 step completed in " << endtime << " seconds." << endl;
    system("pause");
    return 0;
}