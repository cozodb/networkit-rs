
#ifndef NK_CENTRALITY_H
#define NK_CENTRALITY_H

#include "rust/cxx.h"
#include <networkit/centrality/ApproxBetweenness.hpp>
#include <networkit/centrality/ApproxCloseness.hpp>
#include <networkit/centrality/ApproxElectricalCloseness.hpp>

namespace NetworKit
{
    using namespace std;

    inline unique_ptr<ApproxBetweenness> NewApproxBetweenness(
        const Graph &G, double epsilon = 0.01, double delta = 0.1,
        double universalConstant = 1.0)
    {
        return make_unique<ApproxBetweenness>(G, epsilon, delta, universalConstant);
    }
    inline void ApproxBetweennessRanking(ApproxBetweenness &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> ApproxBetweennessScores(ApproxBetweenness &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<ApproxCloseness> NewApproxCloseness(
        const Graph &G, count nSamples, double epsilon, bool normalized,
        uint8_t t)
    {
        ApproxCloseness::ClosenessType type;
        switch (t)
        {
        case 0:
            type = ApproxCloseness::ClosenessType::OUTBOUND;
            break;
        case 1:
            type = ApproxCloseness::ClosenessType::INBOUND;
            break;
        case 2:
            type = ApproxCloseness::ClosenessType::SUM;
            break;
        default:
            type = ApproxCloseness::ClosenessType::OUTBOUND;
        }
        return make_unique<ApproxCloseness>(G, nSamples, epsilon, normalized, type);
    }
    inline void ApproxClosenessRanking(ApproxCloseness &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> ApproxClosenessScores(ApproxCloseness &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<vector<double>> ApproxClosenessGetSquareErrorEstimates(ApproxCloseness &algo)
    {
        return make_unique<vector<double>>(algo.getSquareErrorEstimates());
    }

    inline unique_ptr<ApproxElectricalCloseness> NewApproxElectricalCloseness(
        const Graph &G, double epsilon = 0.1, double kappa = 0.3)
    {
        return make_unique<ApproxElectricalCloseness>(G, epsilon, kappa);
    }
    inline void ApproxElectricalClosenessRanking(ApproxElectricalCloseness &algo, rust::Vec<node> &ks, rust::Vec<double> &vs)
    {
        for (auto &&pair : algo.ranking())
        {
            ks.push_back(pair.first);
            vs.push_back(pair.second);
        }
    }
    inline unique_ptr<vector<double>> ApproxElectricalClosenessScores(ApproxElectricalCloseness &algo)
    {
        return make_unique<vector<double>>(algo.scores());
    }

    inline unique_ptr<vector<double>> ApproxElectricalClosenessComputeExactDiagonal(const ApproxElectricalCloseness &algo, double tol)
    {
        return make_unique<vector<double>>(algo.computeExactDiagonal(tol));
    }

    inline unique_ptr<vector<double>> ApproxElectricalClosenessGetDiagonal(const ApproxElectricalCloseness &algo)
    {
        return make_unique<vector<double>>(algo.getDiagonal());
    }

}

#endif // NK_CENTRALITY_H
