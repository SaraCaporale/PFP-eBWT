#include <stdint.h>
#include <iostream>
#include <set>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include <cmath>
#include <sys/stat.h>
#include <float.h>

#define BUFFERSIZE 1048576

float distSigma(FILE * DA_file, int x, int y) {

  uint64_t sigma=0, x_count=0, y_count=0; 
  unsigned char curr, prev;
  unsigned char* buff = new unsigned char[BUFFERSIZE];

  size_t bytesReadDA=0;
  uint64_t n=0;
  size_t i=0;

	uint64_t c_count=0;

  bytesReadDA=fread(buff, sizeof(unsigned char), BUFFERSIZE, DA_file);
  if (bytesReadDA<=0 && ! feof(DA_file)) {
    std::cerr << "Error reading from file." << std::endl;
    return -1;
  }

  prev=buff[0];

  while (bytesReadDA>=1) {
    for (int i =0; i<bytesReadDA; i++) {
    curr=buff[i];
          
      if (curr==x || curr==y) {
        n++;
        if (curr == prev) {
          c_count++;
          } else {
            if (c_count!=0) {sigma +=(c_count-1);}
            c_count=1; 
            }
            prev=curr;
      }
    }
    bytesReadDA=fread(buff, sizeof(unsigned char), BUFFERSIZE, DA_file);
  }

	sigma+=(c_count-1);
    
  if (!(feof(DA_file))) {
    std::cerr << "Error reading from file." << std::endl;
    return -1;
  }

  return (float)sigma/(n-2);
}

float distRho(FILE * ebwt_file, FILE * DA_file, int x, int y) {
  uint64_t dp=0, nx=0, ny=0;
  uint64_t n=0;

  char prevChar;
  char currChar;
  size_t bytesReadEBWT=0, bytesReadDA=0, i=0, j=0, z=0;
  unsigned char currDA;
  unsigned char* buffDA = new unsigned char[BUFFERSIZE];
  char* buffEBWT = new char[BUFFERSIZE];

  bytesReadEBWT=fread(buffEBWT, sizeof(char), BUFFERSIZE, ebwt_file);
  if (bytesReadEBWT<=0) {
    std::cerr << "Error reading from file." << std::endl;
    return -1;
  }

  prevChar=buffEBWT[0];

  bytesReadDA=fread(buffDA, sizeof(unsigned char), BUFFERSIZE, DA_file);
  if (bytesReadDA<=0) {
    std::cerr << "Error reading from file." << std::endl;
    return -1;
  }

  while(bytesReadEBWT>=1 && bytesReadDA>=1) {
    while (i<bytesReadDA && j<bytesReadEBWT) {
        currDA=buffDA[i];
        currChar=buffEBWT[j];
        if (currDA==x || currDA == y) {
            n++;
            if (currChar == prevChar) {
            if (currDA==x) {nx++;} else {ny++;}
            } else { //the current run ended
            dp+=abs(nx-ny);
            prevChar=currChar;
            if (currDA==x) {nx=1; ny=0;} else {ny=1;nx=0;} 
            }
        }
        i++; j++;
    } 
    if (i==bytesReadDA) {
        bytesReadDA=fread(buffDA, sizeof(unsigned char), BUFFERSIZE, DA_file);
        i=0; 
    } else if (j==bytesReadEBWT){
        bytesReadEBWT=fread(buffEBWT, sizeof(char), BUFFERSIZE, ebwt_file);
        j=0;
    }
  }

  dp+=abs(nx-ny); //last run

  if (ferror(DA_file) || ferror(ebwt_file)) {
    std::cerr << "Error reading from file." << std::endl;
    return -1;
  }
  
  return (float)dp/n;
  
}

float distBWSDm(FILE * DA_file, int x, int y) {
  std::map<uint64_t, uint64_t> m;

  unsigned char prev;
  int j=0;
  unsigned char* buff = new unsigned char[BUFFERSIZE];
  size_t i=0;

  size_t bytesReadDA=0;
  unsigned char currDA;
  while ((bytesReadDA=fread(&currDA, sizeof(unsigned char), 1, DA_file))>=1) { //first run
    if (currDA == x || currDA == y) {
      break;
    }
  }

  if (ferror(DA_file)) {
    std::cerr << "Error reading from file." << std::endl;
    return -1;
  }
  
  prev=currDA;
  uint64_t currRun=1; //first run
  uint64_t n=1;
  uint64_t s=1;

  bytesReadDA=fread(buff, sizeof(unsigned char), BUFFERSIZE, DA_file);
  if (bytesReadDA<=0 && !feof(DA_file)) {
    std::cerr << "Error reading from file." << std::endl;
    return -1;
  }

  while (bytesReadDA>=1) {
    while (i<bytesReadDA) {
      currDA=buff[i];
      if (currDA==x || currDA==y ) {
        n++;
        if (currDA == prev) {
        currRun++;
        } else {
        if (m.find(currRun) != m.end()) {
            m[currRun]++;
        } else {
            m[currRun]=1;
        }
        s++;
        currRun=1;
        }
        prev=currDA;
      }
      i++;
    }
    i=0;
    bytesReadDA=fread(buff, sizeof(unsigned char), BUFFERSIZE, DA_file);
  }

  if (!(feof(DA_file))) {
    std::cerr << "Error reading from file." << std::endl;
    return -1;
  }

  if (m.find(currRun) != m.end()) {
    m[currRun]++;
  } else {
    m[currRun]=1;
  }

  float dm=0;
  for (auto const& pair : m) {
    dm+=((float)pair.first*((float)pair.second/s));
  }

  dm=(dm-1)/((n/2) -1);
  return dm;

}

float distBWSDe(FILE * DA_file, int x, int y) {
  std::map<uint64_t, uint64_t> m;

  unsigned char prev;
  int j=0;
  unsigned char* buff = new unsigned char[BUFFERSIZE];
  size_t i=0;

  size_t bytesReadDA=0;
  unsigned char currDA;
  while ((bytesReadDA=fread(&currDA, sizeof(unsigned char), 1, DA_file))>=1 ) { //first run
    if (currDA == x || currDA == y) {
      break;
    }
  }

  if (ferror(DA_file)) {
    std::cerr << "Error reading from file." << std::endl;
    return -1;
  }

  prev=currDA;
  int currRun=1; //first run
  uint64_t s=1;

  bytesReadDA=fread(buff, sizeof(unsigned char), BUFFERSIZE, DA_file);
  if (bytesReadDA<=0 && !feof(DA_file)) {
    std::cerr << "Error reading from file." << std::endl;
    return -1;
  }

  while (bytesReadDA>=1) {
    while (i<bytesReadDA) {
      currDA=buff[i];
      if (currDA==x || currDA==y ) {
        if (currDA == prev) {
        currRun++;
        } else {
        if (m.find(currRun) != m.end()) {
            m[currRun]++;
        } else {
            m[currRun]=1;
        }
        currRun=1;
        s++;
        }
        prev=currDA;
      }
      i++;
    }
    i=0;
    bytesReadDA=fread(buff, sizeof(unsigned char), BUFFERSIZE, DA_file);
  }

  if (!(feof(DA_file))) {
    std::cerr << "Error reading from file." << std::endl;
    return -1;
  }

  if (m.find(currRun) != m.end()) {
    m[currRun]++;
  } else {
    m[currRun]=1;
  }


  float de=0;
  for (auto const& pair : m) {
    de+=(((float)pair.second/s)*log2((float)pair.second/s));
  }
  
  de=(de*(-1))/log2(s);
  
  return de;

}


int calculate_distances(std::string flag, std::vector<std::string>labels, std::string DA_filename, std::string ebwt_filename) {

  FILE * DA_file;
  if ((DA_file = fopen(DA_filename.c_str(), "r")) == nullptr) {
    std::cerr << "Error: Could not open file '" << DA_filename <<"'." << std::endl;
    perror("open() file failed");
    return -1;
  }

  struct stat file_stat;
  if (stat(DA_filename.c_str(), &file_stat) != 0) {
    std::cerr << "Error getting file stats: " << DA_filename << std::endl;
    return -1;
  }

  FILE * ebwt_file;
  if (flag == "-r") {
    ebwt_file=fopen(ebwt_filename.c_str(), "r");
    if (ebwt_file == nullptr) {
      std::cerr << "Error: Could not open file '" << ebwt_filename <<"'." << std::endl;
      return -1;
    }
  }
  
  float matrix [labels.size()][labels.size()];
  float r;

  for (int i=0; i<labels.size(); i++) {
    matrix [i][i]=0.00;
    for (int j=i+1; j<labels.size(); j++) {
      if (flag == "-s") {
        r=distSigma(DA_file, i, j);    
      } else if (flag == "-r") {
        r=distRho(ebwt_file, DA_file, i, j);
        fseek(ebwt_file, 0, SEEK_SET);
      } else if (flag == "-m") {
        r=distBWSDm(DA_file, i, j); 
      } else if (flag == "-e") {
        r=distBWSDe(DA_file, i, j); 
      }
      
      fseek(DA_file, 0, SEEK_SET);

      matrix[i][j]=r;  
      matrix[j][i]=r;
    }
  }

  fclose(DA_file);


  std::ofstream distMat_file(DA_filename + ".distMat" + (char)toupper(flag.at(1)));
  if (!distMat_file.is_open()) {
    std::cerr << "Error opening file: " << DA_filename << std::endl;
    return {};
  }

  distMat_file << labels.size() << "\n"; 
  for (int i=0; i<labels.size(); i++) {
    distMat_file << labels[i] << " ";
    for (int j=0; j<labels.size(); j++){
      distMat_file << matrix[i][j] << " ";
    }
    distMat_file << "\n";
  }

  distMat_file.close();
  
  return 0;
}

std::vector<std::string> find_labels(std::string labels_filename) { //creates vectors of labels from file with labels after ">"

  std::ifstream labels_file(labels_filename);
  if (!labels_file.is_open()) {
    std::cerr << "Error opening file: " << labels_filename << std::endl;
    return {}; 
  }
  
  std::vector<std::string> labels;
  std::string line;
  int i=0;
  while (std::getline(labels_file, line)) {
    if (line.size() > 0 && line[0] == '>') {
      i++;
      if (line.find(' ') == std::string::npos) {
        line = line.substr(1);
      } else {
        line = line.substr(1, line.find(' ')-1);
      } 
      if (line.size() == 0) {line = "T" + std::to_string(i);}
      labels.push_back(line);
    }
    line.clear();
  }
  
  labels_file.close();

  return labels;
}

int main(int argc, char* argv[]) {

  if (argc < 4) {
    std::cerr << "Usage: " << argv[0] << " [-s | -r | -m | -e] <labels file> <DA filename> <(only with -r) ebwt filename>" << std::endl;
    return -1;
  }

  std::string flag = argv[1];

  if ((flag != "-s") && (flag != "-r") && (flag != "-m") && (flag != "-e")) {
    std::cerr << "Invalid flag: " << flag << std::endl;
    return -1;
  }

  std::string labels_filename = argv[2];
  //Check if second argument is a file
  struct stat file_stat_labels;
  if (stat(labels_filename.c_str(), &file_stat_labels) != 0) {
    std::cerr << "Error: File '" << labels_filename << "' does not exist." << std::endl;
    return -1;
  }
  std::vector<std::string> labels = find_labels(labels_filename);

  std::string DA_filename = argv[3];
  //Check if third argument is a file
  struct stat file_stat_DA;
  if (stat(DA_filename.c_str(), &file_stat_DA) != 0) {
    std::cerr << "Error: File '" << DA_filename << "' does not exist." << std::endl;
    return -1;
  }

  std::string ebwt_filename;
  //ebwt file needed only with -r flag
  if (flag=="-r") {
    if (argc < 5) {
      std::cerr << "Error: to calculate rho distance you need to input ebwt file too." << std::endl;
      return -1;
    }
    ebwt_filename = argv[4];
    struct stat file_stat_ebwt; 
    if (stat(ebwt_filename.c_str(), &file_stat_ebwt) != 0) {
      std::cerr << "Error: File '" << ebwt_filename << "' does not exist." << std::endl;
      return -1;
    }
  }


  int r = calculate_distances(flag, labels, DA_filename, ebwt_filename);
  if (r==0) {
    std::cout << "Distances computed successfully." << std::endl;
    return 0;
  } else {
    std::cerr << "Error: Failed to compute distances." << std::endl;
    return 1;
  }
  
}
