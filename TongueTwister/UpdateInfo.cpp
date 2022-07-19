#include <iomanip>
#include <iostream>
#include "UpdateInfo.hpp"


UpdateInfo::UpdateInfo(void) {

}

UpdateInfo::~UpdateInfo(void) {

    for (std::map<std::string,AcceptTries*>::iterator it = info.begin(); it != info.end(); it++)
        delete it->second;
    info.clear();
}

void UpdateInfo::accept(std::string key) {

    AcceptTries* at = getUpdateInfo(key);
    at->numTries++;
    at->numAccepts++;
}

AcceptTries* UpdateInfo::getUpdateInfo(std::string key) {

    std::map<std::string,AcceptTries*>::iterator it = info.find(key);
    if (it == info.end())
        {
        AcceptTries* newInfo = new AcceptTries;
        info.insert( std::make_pair(key,newInfo) );
        return newInfo;
        }
    return it->second;
}

void UpdateInfo::print(void) {

    int len = 0;
    for (std::map<std::string,AcceptTries*>::iterator it = info.begin(); it != info.end(); it++)
        {
        if (it->first.length() > len)
            len = (int)it->first.length();
        }
        
    for (std::map<std::string,AcceptTries*>::iterator it = info.begin(); it != info.end(); it++)
        {
        std::cout << "   * Acceptance rate for update of " << it->first << " ";
        for (int i=0; i<len-it->first.length(); i++)
            std::cout << " ";
        std::cout << "= ";
        if (it->second->numTries > 0)
            std::cout << std::fixed << std::setprecision(1) << ((double)it->second->numAccepts / it->second->numTries) * 100.0 << "%";
        else
            std::cout << "N/A";

        std::cout << std::endl;
        }
}

void UpdateInfo::reject(std::string key) {

    AcceptTries* at = getUpdateInfo(key);
    at->numTries++;
}

