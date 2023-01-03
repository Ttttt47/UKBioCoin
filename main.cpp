#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <limits>
#include <Eigen/Dense>
#include <Eigen/Core>

using std::vector;
using namespace std;
 
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

void Read_matrix_table(string filename, vector<string> colnames=NULL)
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
    getline(fin, line);  // skip first line.
    while(getline(fin, line))
    {   
        istringstream iline(line);
        iline >> word;  // ignore first colunm.
        while (iline >> word)
        {   
            if (word == "NA")
            {
                el = std::numeric_limits<double>::quiet_NaN();
            }
            else
            {
                el = std::stod(word);
            }
            line_data.push_back(el);
        }
        for (int i = 0; i < line_data.size(); ++i) 
        {
            std::cout << line_data[i] << " ";
        }
        cout << endl;
        line_data = vector<double>();
    }
    exit(0);
}


int main()
{   
    string cov_xy_filename = "./test_data/cov_xy.table";
    string var_x_filename = "./test_data/var_x_sam.table";
    string cov_yy_filename = "./test_data/cov_yy.table";
    vector<string> covar = {"X31.0.0", "X34.0.0", "X21000.0.0", 
                            "PC1", "PC2", "PC3", "PC4", "PC5"};
    Read_matrix_table(cov_yy_filename);
    return 0;
}