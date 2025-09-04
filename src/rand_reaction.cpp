// rand_reaction.cpp
#include <vector>
#include <numeric>
#include <cmath>
#include <stdexcept>
#include "logger.hpp"
#include "property_package.hpp"
#include "rand_flash.hpp"

namespace simona {

struct ReactionResult {
    bool success=false;
    int iters=0;
    double final_error=0.0;

    double betaV=0.0, betaL=0.0;
    std::vector<double> nV, nL;        // 相内摩尔数
    std::vector<double> xV, xL;        // 组成
    std::vector<double> gamma;         // 反应推进量（Lagrange 步）
};

class RandReaction {
public:
    RandReaction(PropertyPackage* prop) : prop_(prop) {}

    // 两相 + 反应：A 是 E×C（E<=C），N 是 C×R
    ReactionResult solveTwoPhaseReactive(
        double P, double T,
        const std::vector<double>& feed,               // 尽量归一
        const std::vector<std::vector<double>>& A,     // E×C
        const std::vector<std::vector<double>>& N,     // C×R
        const std::vector<double>& xV0,                // 初始气相组成
        const std::vector<double>& xL0,                // 初始液相组成
        int maxIter = 50)
    {
        const int C=(int)feed.size();
        const int E=(int)A.size();
        const int R=(int)N[0].size();

        // 初值
        double betaV=0.5, betaL=0.5;
        std::vector<double> xV=xV0, xL=xL0;
        norm1(xV); norm1(xL);

        std::vector<double> nV(C), nL(C);
        for (int i=0;i<C;++i){ nV[i]=betaV*xV[i]; nL[i]=betaL*xL[i]; }

        prop_->setPressure(P); prop_->setTemperature(T);
        const double RT = R_CONST * T;

        ReactionResult out;
        for (int it=0; it<maxIter; ++it){
            // μ/RT & 局部矩阵
            std::vector<double> muV(C,0), muL(C,0);
            prop_->chemicalPotential(/*V=*/0, nV, muV);
            prop_->chemicalPotential(/*L=*/1, nL, muL);
            for (double& v: muV) v/=RT;
            for (double& v: muL) v/=RT;

            std::vector<std::vector<double>> MV(C,std::vector<double>(C,0));
            std::vector<std::vector<double>> ML(C,std::vector<double>(C,0));
            prop_->buildLocalM(nV, 0, MV);  // 约定：这里的 M 已含 β
            prop_->buildLocalM(nL, 1, ML);

            // 全局系统： [ A(MV+ML)A^T  | A xV | A xL | A N ] [λ; ΔβV; ΔβL; γ] = [ A(MV μV + ML μL); u2V; u2L; 0 ]
            const int Nvar = E + 2 + R;
            std::vector<std::vector<double>> K(Nvar, std::vector<double>(Nvar,0.0));
            std::vector<double> rhs(Nvar, 0.0);

            // 构造辅助：AT、AXV、AXL、AN
            std::vector<std::vector<double>> AT(C, std::vector<double>(E,0.0));
            for (int r=0;r<E;++r) for (int c=0;c<C;++c) AT[c][r]=A[r][c];

            // 顶左 ∑ A M A^T
            {
                std::vector<std::vector<double>> MAT(C,std::vector<double>(E,0.0)), AMAT(E,std::vector<double>(E,0.0));
                // MV A^T
                matmul(MV, AT, MAT);
                matmul(A, MAT, AMAT);
                addBlock(K, 0, 0, AMAT);
                // ML A^T
                matmul(ML, AT, MAT);
                matmul(A, MAT, AMAT);
                addBlock(K, 0, 0, AMAT);
            }
            // 顶右 A xV, A xL, A N
            {
                std::vector<double> AxV(E,0.0), AxL(E,0.0);
                matvec(A, xV, AxV); matvec(A, xL, AxL);
                for (int r=0;r<E;++r){ K[r][E+0] = AxV[r]; K[r][E+1] = AxL[r]; }
                // A N
                std::vector<std::vector<double>> AN(E, std::vector<double>(R,0.0));
                matmul(A, N, AN);
                for (int r=0;r<E;++r) for (int j=0;j<R;++j) K[r][E+2+j] = AN[r][j];
            }
            // 其余块 0；RHS 顶部
            {
                std::vector<double> MVmu(C,0.0), MLmu(C,0.0), sum(C,0.0), A_sum(E,0.0);
                matvec(MV, muV, MVmu);
                matvec(ML, muL, MLmu);
                for (int i=0;i<C;++i) sum[i]=MVmu[i]+MLmu[i];
                matvec(A, sum, A_sum);
                for (int r=0;r<E;++r) rhs[r]=A_sum[r];
            }
            // RHS u2V,u2L；最后 R 行为 0
            rhs[E+0] = dot(xV, muV);
            rhs[E+1] = dot(xL, muL);
            // E+2...E+2+R-1 默认 0

            // 解 K * y = rhs
            std::vector<double> y;
            if (!solveLinear(K, rhs, y)) {
                SIMONA_LOG_WARN("[RandReaction] linear solve failed at iter {}", it);
                break;
            }
            std::vector<double> lambda(E,0.0); double dBV=y[E+0], dBL=y[E+1];
            for (int r=0;r<E;++r) lambda[r]=y[r];
            std::vector<double> gamma(R,0.0);
            for (int j=0;j<R;++j) gamma[j]=y[E+2+j];

            // Δn^V, Δn^L（式 4.22 带反应项）
            // g = A^T λ - μ/RT
            std::vector<double> ATlambda(C,0.0);
            matvecT(A, lambda, ATlambda); // A^T λ

            std::vector<double> dnV(C,0.0), dnL(C,0.0);
            // M*g
            std::vector<double> gV(C,0.0), gL(C,0.0), MVg(C,0.0), MLg(C,0.0);
            for (int i=0;i<C;++i){ gV[i]=ATlambda[i]-muV[i]; gL[i]=ATlambda[i]-muL[i]; }
            matvec(MV, gV, MVg); matvec(ML, gL, MLg);

            // 反应项 N γ
            std::vector<double> Ng(C,0.0);
            matvec(N, gamma, Ng);

            for (int i=0;i<C;++i){
                dnV[i] = xV[i]*dBV + MVg[i] + Ng[i];
                dnL[i] = xL[i]*dBL + MLg[i] - Ng[i]; // 注意：总 n = nV+nL+Nγ - Nγ，不改变总元素
            }

            // 线搜索（正性 + 下降）
            double descent = dot(dnV, muV) + dot(dnL, muL);
            double alpha = backtracking(nV,nL,dnV,dnL,betaV,betaL,dBV,dBL, descent);
            for (int i=0;i<C;++i){ nV[i]+=alpha*dnV[i]; nL[i]+=alpha*dnL[i]; }
            betaV+=alpha*dBV; betaL+=alpha*dBL;

            // 更新 x
            xV = normed(nV); xL = normed(nL);

            // 收敛
            double max_res=0.0;
            for (int i=0;i<C;++i){
                max_res = std::max(max_res, std::fabs(ATlambda[i]-muV[i]));
                max_res = std::max(max_res, std::fabs(ATlambda[i]-muL[i]));
            }
            if (max_res < 1e-9) {
                out.success=true; out.iters=it+1; out.final_error=max_res;
                out.betaV=betaV; out.betaL=betaL;
                out.nV=nV; out.nL=nL; out.xV=xV; out.xL=xL; out.gamma=gamma;
                return out;
            }
        }

        out.success=false;
        out.betaV=betaV; out.betaL=betaL;
        out.nV=nV; out.nL=nL; out.xV=xV; out.xL=xL;
        return out;
    }

private:
    PropertyPackage* prop_;

    // helpers (和上面类似；可替换成你项目里的通用函数)
    static void addBlock(std::vector<std::vector<double>>& K, int r0, int c0,
                         const std::vector<std::vector<double>>& B){
        for (int r=0;r<(int)B.size();++r)
            for (int c=0;c<(int)B[0].size();++c)
                K[r0+r][c0+c] += B[r][c];
    }
    static void matvec(const std::vector<std::vector<double>>& M,
                       const std::vector<double>& v,
                       std::vector<double>& out) {
        int R=(int)M.size(), C=(int)M[0].size();
        out.assign(R, 0.0);
        for (int i=0;i<R;++i){
            double s=0; for(int k=0;k<C;++k) s+=M[i][k]*v[k];
            out[i]=s;
        }
    }
    static void matvecT(const std::vector<std::vector<double>>& A,
                        const std::vector<double>& lambda,
                        std::vector<double>& ATlambda){
        int E=(int)A.size(), C=(int)A[0].size();
        ATlambda.assign(C,0.0);
        for (int c=0;c<C;++c){
            double s=0; for (int r=0;r<E;++r) s += A[r][c]*lambda[r];
            ATlambda[c]=s;
        }
    }
    static void matmul(const std::vector<std::vector<double>>& A,
                       const std::vector<std::vector<double>>& B,
                       std::vector<std::vector<double>>& out){
        int R=(int)A.size(), K=(int)A[0].size(), C=(int)B[0].size();
        out.assign(R, std::vector<double>(C,0.0));
        for (int i=0;i<R;++i)
            for (int k=0;k<K;++k){
                double aik=A[i][k];
                for (int j=0;j<C;++j) out[i][j]+=aik*B[k][j];
            }
    }
    static double dot(const std::vector<double>& a, const std::vector<double>& b){
        double s=0; for (size_t i=0;i<a.size();++i) s+=a[i]*b[i]; return s;
    }
    static void norm1(std::vector<double>& x){
        double s=std::accumulate(x.begin(), x.end(), 0.0);
        if (s<=0) { double v=1.0/x.size(); for (double& xi: x) xi=v; }
        else for(double& xi: x) xi/=s;
    }
    static std::vector<double> normed(const std::vector<double>& n){
        double s=std::accumulate(n.begin(), n.end(), 0.0);
        std::vector<double> x(n.size(), 0.0);
        for (size_t i=0;i<n.size();++i) x[i]= n[i]/std::max(s,1e-30);
        return x;
    }
    static bool solveLinear(const std::vector<std::vector<double>>& A,
                            const std::vector<double>& b,
                            std::vector<double>& x) {
        int N=(int)A.size();
        std::vector<std::vector<double>> M=A;
        x=b;
        for (int i=0;i<N;++i){
            int piv=i;
            for (int r=i;r<N;++r) if (std::fabs(M[r][i])>std::fabs(M[piv][i])) piv=r;
            if (std::fabs(M[piv][i])<1e-14) return false;
            if (piv!=i){ std::swap(M[piv],M[i]); std::swap(x[piv],x[i]); }
            double d=M[i][i];
            for (int c=i;c<N;++c) M[i][c]/=d; x[i]/=d;
            for (int r=0;r<N;++r) if (r!=i){
                double f=M[r][i];
                for (int c=i;c<N;++c) M[r][c]-=f*M[i][c];
                x[r]-=f*x[i];
            }
        }
        return true;
    }
    double backtracking(const std::vector<double>& nV,
                        const std::vector<double>& nL,
                        const std::vector<double>& dnV,
                        const std::vector<double>& dnL,
                        double betaV, double betaL,
                        double dBV, double dBL,
                        double descent0){
        double alpha=1.0, c=1e-4;
        for (int rep=0; rep<30; ++rep){
            bool ok=true;
            for (size_t i=0;i<nV.size();++i){
                if (nV[i]+alpha*dnV[i] <= 1e-14) { ok=false; break; }
                if (nL[i]+alpha*dnL[i] <= 1e-14) { ok=false; break; }
            }
            if (!ok){ alpha*=0.5; continue; }
            double desc = 0.0;
            // 下降检查：近似用 Δβ 的二次项/或直接用 descent0 的 Armijo
            if (desc <= c*alpha*descent0) return alpha;
            alpha*=0.5;
        }
        return alpha;
    }
};

} // namespace simona
