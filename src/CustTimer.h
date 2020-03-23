#ifndef CUSTTIMER_H
#define CUSTTIMER_H

/**

This class presents a basic timer that is used in KLUSolver to know on which part of the solver
most time were taken.

**/
class CustTimer{
    public:
        CustTimer():start_(std::chrono::system_clock::now()){
            end_ = start_;
        };

        double duration(){
            end_ = std::chrono::system_clock::now();
            std::chrono::duration<double> res = end_ - start_;
            return res.count();
        }
    private:
        std::chrono::time_point<std::chrono::system_clock> start_;
        std::chrono::time_point<std::chrono::system_clock> end_;
};

#endif //CUSTTIMER_H