# include "precise.hpp"



precise_real_type operator ""_p(const char* a){
	return static_cast<precise_real_type>(a);
}

// //Old implementation;
// precise_real_type operator ""_p(const char* a){
// 	precise_real_type res=0;
// 	std::string number=a;

// 	std::size_t found2=number.find("e");
// 	precise_real_type fac=1;

// 	if(found2!=std::string::npos){
// 		fac=precise_real_type(std::stoi(number.substr(found2+1,number.size()-found2-1)));
// 		fac=pow(precise_real_type(10),fac);
// 		number=number.substr(0,found2);
// 	}

// 	std::size_t found=number.find(".");

// 	if(found!=std::string::npos && found!=number.size()-1){
// 		for(int i=number.substr(0,found).size()-1; i>=0; i--){
// 			res+=precise_real_type(std::stoi(number.substr(i,1)))*pow(precise_real_type(10), number.substr(0,found).size()-1-i);
// 		}
// 		for(int i=0; i<=number.substr(found+1).size()-1;i++){
// 			res+=precise_real_type(std::stoi(number.substr(found+1+i,1)))*pow(precise_real_type(10),-i-1);
// 		}
// 	}
// 	else if (found==number.size()-1){
// 		number=number.substr(0,number.size()-1);
// 		for(int i=number.size()-1; i>=0; i--){
// 			res+=precise_real_type(std::stoi(number.substr(i,1)))*pow(precise_real_type(10), number.size()-1-i);
// 		}
// 	}
// 	else{
// 		for(int i=number.size()-1; i>=0; i--){
// 			res+=precise_real_type(std::stoi(number.substr(i,1)))*pow(precise_real_type(10), number.size()-1-i);
// 		}
// 	}
// 	return res*fac;
// }