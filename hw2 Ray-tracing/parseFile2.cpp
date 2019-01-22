//File parsing example

#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(){
  string line;

  string fileName = "spheres1.scn";

  // open the file containing the scene description
  ifstream input(fileName);

  // check for errors in opening the file
  if(input.fail()){
    cout << "Can't open file '" << fileName << "'" << endl;
    return 0;
  }
  
  // determine the file size (this is optional -- feel free to delete the 6 lines below)
  streampos begin,end;
  begin = input.tellg();
  input.seekg(0, ios::end);
  end = input.tellg();
  cout << "File '" << fileName << "' is: " << (end-begin) << " bytes long.\n\n";
  input.seekg(0, ios::beg);

  
  //Loop through reading each line
  string command;
  while(input >> command) { //Read first word in the line (i.e., the command type)
    
    if (command[0] == '#'){
      getline(input, line); //skip rest of line
      cout << "Skipping comment: " << command  << line <<  endl;
      continue;
    }
    
    
    if (command == "sphere"){ //If the command is a sphere command
       float x,y,z,r;
       input >> x >> y >> z >> r;
       printf("Sphere as position (%f,%f,%f) with radius %f\n",x,y,z,r);
    }
    else if (command == "background"){ //If the command is a background command
       float r,g,b;
       input >> r >> g >> b;
       printf("Background color of (%f,%f,%f)\n",r,g,b);
    }
    else if (command == "output_image"){ //If the command is an output_image command
       string outFile;
       input >> outFile;
       printf("Render to file named: %s\n", outFile.c_str());
    }
    else {
      getline(input, line); //skip rest of line
      cout << "WARNING. Do not know command: " << command << endl;
    }
  }
  
  return 0;
}
