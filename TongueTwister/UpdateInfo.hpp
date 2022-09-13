#ifndef UpdateInfo_hpp
#define UpdateInfo_hpp

#include <map>
#include <string>



struct AcceptTries {

            AcceptTries(void) { numTries = 0; numAccepts = 0; }
    int     numTries;
    int     numAccepts;
};

class UpdateInfo {

    public:
        static UpdateInfo&                  updateInfo(void) {
                                                static UpdateInfo uiPool;
                                                return uiPool;
                                            }
        void                                accept(std::string& key);
        void                                print(void);
        void                                reject(std::string& key);
        AcceptTries*                        getUpdateInfo(std::string& key);
        
    private:
                                            UpdateInfo(void);
                                           ~UpdateInfo(void);
        std::map<std::string,AcceptTries*>  info;
};

#endif
