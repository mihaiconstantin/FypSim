#include <Rcpp.h>
#include "Statistics.h"


// Simulate data for the GRM based on a vector of
// abilities and a matrix of item parameters.
Rcpp::NumericMatrix Statistics::SimulateGrm(const Rcpp::DoubleVector &theta, const Rcpp::NumericMatrix &item_parameters)
{
    // Get some idea about the size of the arguments passed in.
    int total_rows = (int) theta.length();
    int total_cols = item_parameters.nrow();
    int total_steps = item_parameters.ncol() - 1;


    // Let's build what we are after. We'll populate this accordingly.
    Rcpp::NumericMatrix item_means(total_rows, total_cols);
    Rcpp::NumericMatrix item_scores(total_rows, total_cols);


    // First looping over the items, then looping over the steps.
    // We'll relay on some vector operations to gain speed.
    for (int item = 0; item < total_cols; ++item)
    {
        // Let's define the structures that will hold data within
        // the loop that goes over each of the items steps.
        Rcpp::NumericMatrix item_step_response_functions(total_rows, total_steps);
        Rcpp::NumericMatrix scores(total_rows, total_steps);


        // Also, we need to take a number of random draws that we'll use
        // to compare the isrf to in order to generate the response
        // values. The number must be equal to the length of the
        // ability estimate.
        Rcpp::DoubleVector random_draws = Rcpp::runif(total_rows, 0, 1);


        // Starting to iterate over each step.
        for (int step = 0; step < total_steps; ++step)
        {
            Rcpp::DoubleVector logit = item_parameters(item, 0) * (theta - item_parameters(item, step + 1));
            Rcpp::DoubleVector isrf = 1 / (1 + Rcpp::exp(-logit));


            // Now that we've completed the computations, let's store the
            // results for this step in the item matrices allocated to
            // the loop that goes over the item steps.
            item_step_response_functions(Rcpp::_, step) = isrf;
            scores(Rcpp::_, step) = (isrf - random_draws) >= 0;
        }


        // Before we move to the next item, we must sum all the columns
        // for our item step response functions and scores matrices.
        // This will result in the expected item means for the
        // current item and the item values, respectively.
        // We store everything into the matrices built
        // for the loop that goes over the item.
        item_means(Rcpp::_, item) = Rcpp::rowSums(item_step_response_functions);
        item_scores(Rcpp::_, item) = Rcpp::rowSums(scores);
    }

    // Let's return a nice array of numeric matrices, where
    // the first index holds the expected item means and
    // the second index the data points.
    Rcpp::NumericMatrix result[2] = {item_means, item_scores};

    return result[1];
}



// Call the estimateParametersGrm R function that serves
// as a wrapper for the mirt package. At the end make
// sure your return the matrix of item parameters.
// For a detailed description on the output,
// check the documentation within file
// R/ParametersEstimation.R
Rcpp::NumericMatrix Statistics::EstimateParametersGrm(const Rcpp::NumericMatrix &data)
{
    Rcpp::Environment FypSim_env("package:FypSim");
    Rcpp::Function estimateParametersGrm = FypSim_env["estimateParametersGrm"];
    return estimateParametersGrm(data);
}



// region Theta WML and PolyLzStar (response vector)

// Wilco's device for WLE theta estimation. Runs for each
// item response vector at a time.
Rcpp::NumericVector Statistics::EstimateThetaWml(const Rcpp::NumericVector &data, const Rcpp::NumericMatrix &ipar, double thStart, double tolerance)
{
    // TODO: Why do you initialize them to 10?
    // Declarations.
    Rcpp::NumericVector thHat(4);
    Rcpp::NumericVector P(10), Q(10);

    double a, b;
    double delta = 1;
    double Pr;
    double logit, LogLD0 = 0, LogLD1 = 0, LogLD2 = 0;
    int m, j, s, t = 0;
    int converg = 1;
    int valOut = 0;

    double D0, D1, D2, D3;
    double info, d1info, jacob, d1jacob;
    double scoreFunc, D1scoreFunc;

    double t1, t2, t11, t12, t21, t22;
    double P1ac, P2ac, Q1ac, Q2ac;


    // Counters.
    int K = ipar.ncol();
    int M = K - 1;
    int J = ipar.nrow();
    int MaxIter = 100;


    // Initialize.
    double thEst = thStart;


    // Newton-Raphson method for root approximation.
    while ((abs(delta) > tolerance) & (t < MaxIter))
    {
        t = t + 1;

        LogLD1 = 0;
        LogLD2 = 0;
        LogLD0 = 0;
        info = d1info = jacob = d1jacob = 0;


        // Looping over the items.
        for (j = 0; j < J; j++) {

            // Step probabilities.
            a = ipar(j, 0);
            P[0] = 1;
            Q[0] = 0;

            for (m = 0; m < M; m++)
            {
                b = ipar(j, (m + 1));
                logit = a * (thEst - b);
                P[m + 1] = exp(logit) / (1 + exp(logit));
                Q[m + 1] = 1 - P[m + 1];
            }

            P[M + 1] = 0;
            Q[M + 1] = 1;


            // Log-Likelihood and derivatives.
            s = data[j];
            Pr = P[s] - P[(s + 1)];
            LogLD0 = LogLD0 + log(Pr);
            LogLD1 = LogLD1 + a * (1 - P[s] - P[(s + 1)]);
            LogLD2 = LogLD2 - pow(a, 2) * (P[s] * (1 - P[s]) + P[(s + 1)] * (1 - P[(s + 1)]));


            // Compute Derivatives and Information Function values, and Likelihood:
            for (m = 0; m < (M + 1); m++)
            {
                // Orc probability.
                D0 = P[m] - P[m + 1];

                // First order derivatives.
                D1 = a * P[m] * Q[m] - a * P[m + 1] * Q[m + 1];

                // Second order derivatives.
                t1 = pow(a, 2) * P[m] * Q[m] * (Q[m] - P[m]);
                t2 = pow(a, 2) * P[m + 1] * Q[m + 1] * (Q[m + 1] - P[m + 1]);
                D2 = t1 - t2;

                // Third order derivatives.
                P1ac = a * P[m] * Q[m];
                Q1ac = -a * P[m] * Q[m];
                P2ac = a * P[m + 1] * Q[m + 1];
                Q2ac = -a * P[m + 1] * Q[m + 1];
                t11 = pow(a, 2) * P1ac * pow(Q[m], 2) + pow(a, 2) * P[m] * 2 * Q[m] * Q1ac;
                t12 = pow(a, 2) * 2 * P[m] * P1ac * Q[m] + pow(a, 2) * pow(P[m], 2) * Q1ac;
                t21 = pow(a, 2) * P2ac * pow(Q[m + 1], 2) + pow(a, 2) * P[m + 1] * 2 * Q[m + 1] * Q2ac;
                t22 = pow(a, 2) * 2 * P[m + 1] * P2ac * Q[m + 1] + pow(a, 2) * pow(P[m + 1], 2) * Q2ac;
                D3 = (t11 - t12) - (t21 - t22);

                // Test.
                info = info + pow(D1, 2) / D0;
                d1info = d1info + D1 * ((2 * D2 / D0) - (pow(D1, 2) / pow(D0, 2)));
                jacob = jacob + D1 * D2 / D0;
                t1 = (pow(D2, 2) + (D1 * D3)) * D0;
                t2 = pow(D1, 2) * D2;
                d1jacob = d1jacob + (t1 - t2) / pow(D0, 2);

            } // End of Derivatives and Information Function values, and Likelihood.

        } // End looping over the items.


        // Update estimate WML.
        scoreFunc = jacob / (2 * info) + LogLD1;
        t1 = info - info * d1jacob - pow(jacob, 2);
        t2 = pow(2 * info, 2);
        D1scoreFunc = t1 / t2 + LogLD2;
        delta = (scoreFunc / D1scoreFunc);
        thEst = thEst - delta;

    }  // End of Newton-Raphson method.


    // Compute final stats.
    info = 0;

    // Looping over the items.
    for (j = 0; j < J; j++)
    {
        // Step probabilities.
        a = ipar(j, 0);
        P[0] = 1;
        Q[0] = 0;

        for (m = 0; m < M; m++)
        {
            b = ipar(j, (m + 1));
            logit = a * (thEst - b);
            P[m + 1] = exp(logit) / (1 + exp(logit));
            Q[m + 1] = 1 - P[m + 1];
        }

        P[M + 1] = 0;
        Q[M + 1] = 1;

        for (m = 0; m < (M + 1); m++)
        {
            // Orc probability.
            D0 = P[m] - P[m + 1];

            // First order derivatives.
            D1 = a * P[m] * Q[m] - a * P[m + 1] * Q[m + 1];
            info = info + pow(D1, 2) / D0;
        }

    } // End looping over the items.


    // Prepare the return.
    if (t < MaxIter) { converg = 1; } else { converg = 0; }

    if(abs((int) thEst) < 4 && thEst == thEst) { valOut = 1; } else { valOut = 0; }


    thHat[0] = thEst;
    thHat[1] = 1 / sqrt(info);
    thHat[2] = converg;
    thHat[3] = valOut;


    return (thHat);
}



// Wilco's device for polychotomous corrected LZ estimation.
// Runs for each item response vector at a time.
Rcpp::NumericVector Statistics::EstimatePolyLzStar(const Rcpp::NumericVector &data, const Rcpp::NumericMatrix &ipar, const double &theta)
{
    // TODO: The initialization at 10 means that envision a maximum number of 10 item steps.
    // Declarations.
    Rcpp::NumericVector lz(2), P(10), Q(10);

    lz[0] = theta;
    const int maxJ = 80;
    const int maxK = 10;

    int K = ipar.ncol();
    int M = K - 1;
    int J = ipar.nrow();
    int m, j, l;
    int isc;

    double a, b, logit;
    double info, jacob;

    Rcpp::NumericMatrix D0(maxJ, maxK), D1(maxJ, maxK), D2(maxJ, maxK), sij(maxJ, maxK), wtilde(maxJ, maxK);
    Rcpp::NumericMatrix Dmat(maxK, maxK);

    double t1, t2;
    double s0, cntel, cnnoem, cn;
    double varLz;

    Rcpp::NumericVector varLzT1(maxK);

    double stat;

    // Compute ORCs and Derivatives.
    cntel = cnnoem = info = jacob = 0;

    // Looping over items.
    for (j = 0; j < J; j++)
    {
        a = ipar(j, 0);
        P[0] = 1;
        Q[0] = 0;

        for (m = 0; m < M; m++)
        {
            b = ipar(j, (m + 1));
            logit = a * (theta - b);
            P[m + 1] = exp(logit) / (1 + exp(logit));
            Q[m + 1] = 1 - P[m + 1];
        }

        P[M + 1] = 0;
        Q[M + 1] = 1;

        for (m = 0; m < K; m++)
        {
            D0(j, m) = P[m] - P[m + 1];

            // First order derivatives.
            D1(j, m) = a * (P[m] * Q[m] - P[m + 1] * Q[m + 1]);

            // Second order derivatives.
            t1 = pow(a, 2) * P[m] * Q[m] * (Q[m] - P[m]);
            t2 = pow(a, 2) * P[m + 1] * Q[m + 1] * (Q[m + 1] - P[m + 1]);
            D2(j, m) = t1 - t2;

            info = info + pow(D1(j, m), 2) / D0(j, m);
            jacob = jacob + D1(j, m) * D2(j, m) / D0(j, m);
            sij(j, m) = D1(j, m) / D0(j, m);
            cntel = cntel + D1(j, m) * log(D0(j, m));
            cnnoem = cnnoem + D1(j, m) * sij(j, m);
        }

    } // End looping over items.


    s0 = jacob / (2 * info);
    cn = cntel / cnnoem;


    for (j = 0; j < J; j++)
    {
        for (m = 0; m < K; m++)
        {
            wtilde(j, m) = log(D0(j, m)) - cn * sij(j, m);
        }
    }


    varLz = 0;
    stat = 0;

    // Lopping over items.
    for (j = 0; j < J; j++)
    {
        for (m = 0; m < K; m++)
        {
            if (data[j] == m) { isc = 1; } else { isc = 0; }
            stat = stat + (isc - D0(j, m)) * log(D0(j, m));
        }

        for (m = 0; m < K; m++)
        {
            for (l = 0; l < K; l++)
            {
                Dmat(m, l) = -D0(j, m) * D0(j, l);
            }
        }

        for (m = 0; m < K; m++)
        {
            Dmat(m, m) = D0(j, m) * (1 - D0(j, m));
        }

        for (m = 0; m < K; m++)
        {
            varLzT1[m] = 0;
        }

        for (m = 0; m < K; m++)
        {
            for (l = 0; l < K; l++)
            {
                varLzT1[m] = varLzT1[m] + wtilde(j, l) * Dmat(l, m);
            }
        }

        for (m = 0; m < K; m++)
        {
            varLz = varLz + varLzT1[m] * wtilde(j, m);
        }

    } // End loop across items.


    // Prepare the return.
    lz[1] = (stat + cn * s0) / sqrt(varLz);


    // region TODO: debug scenarios when we use more than 70 items
    //std::cout << "varLz: " << varLz << std::endl;
    //std::cout << "stat: " << stat << std::endl;
    //std::cout << "cn: " << cn << std::endl;
    //std::cout << "s0: " << s0 << std::endl;
    //std::cout << "lz: " << (stat + cn * s0) / sqrt(varLz) << std::endl;
    //std::cout << "vector: " << data << std::endl;
    //std::cout << "theta: " << theta << std::endl;
    //std::cout << "------------------" << std::endl;
    // endregion


    return (lz);
}



// Wilco's device in case the WML estimated with the
// Newton-Raphson optimization method does not
// converge.
Rcpp::NumericVector Statistics::EstimateThetaBisection(Rcpp::NumericVector data, Rcpp::NumericMatrix ipar, double estmeth, double thStart, double tolerance, double bislowin, double bisupin, int nbisin)
{
    // Declarations.
    Rcpp::NumericVector thHat(3);

    int m, j, nbis, ep;
    double a, b, logit;
    double th, tha, thb, thc;
    double Pr, w1, w2, LogLD1;
    double InfoVal, Jacob;
    double D1, D2;
    double t1, t2;

    int s;
    int J = ipar.nrow();
    int K = ipar.ncol();
    int M = K - 1;

    Rcpp::NumericVector P(10), Q(10), funcval(3);

    tha = bislowin, thb = bisupin;
    thc = 0;


    // Start the bisections.
    for (nbis = 0; nbis < nbisin; nbis++)
    {
        // Evaluation points loop.
        for (ep = 0; ep < 3; ep++)
        {
            // Choose point.
            if (ep == 0) { th = tha; }
            if (ep == 1) { th = thb; }
            if (ep == 2) { th = thc; }


            // Compute first derivative and log-likelihood.
            LogLD1 = 0;
            InfoVal = 0;
            Jacob = 0;


            // Loop across items.
            for (j = 0; j < J; j++)
            {
                // Step probabilities.
                a = ipar(j, 0);
                P[0] = 1;
                Q[0] = 0;

                for (m = 0; m < M; m++)
                {
                    b = ipar(j, (m + 1));
                    logit = a * (th - b);
                    P[m + 1] = exp(logit) / (1 + exp(logit));
                    Q[m + 1] = 1 - P[m + 1];
                }

                P[M + 1] = 0;
                Q[m + 1] = 1;

                // Derivative.
                s = data[j];
                Pr = P[s] - P[s + 1];
                w1 = P[s] * (1 - P[s]);
                w2 = P[s + 1] * (1 - P[s + 1]);
                LogLD1 = LogLD1 + a * ((w1 - w2) / Pr);


                // Information values and Jacobian method.
                for (m = 0; m < (M + 1); m++)
                {
                    // First order derivatives.
                    Pr = P[m] - P[(m + 1)];
                    D1 = a * (P[m] * (1 - P[m]) - P[m + 1] * (1 - P[m + 1]));
                    InfoVal = InfoVal + pow(D1, 2) / Pr;

                    // Second order derivatives.
                    t1 = pow(a, 2) * P[m] * Q[m] * (Q[m] - P[m]);
                    t2 = pow(a, 2) * P[m + 1] * Q[m + 1] * (Q[m + 1] - P[m + 1]);
                    D2 = t1 - t2;
                    Jacob = Jacob + D1 * D2 / Pr;
                }

            } // End loop across items.

            funcval[ep] = LogLD1 + estmeth * Jacob / (2 * InfoVal);

        } // End evaluation points.


        if ((funcval[0] > 0) & (funcval[2] <= 0))
        {
            thb = thc;
            thc = (tha + thb) / 2;
        }


        if ((funcval[2] >= 0) & (funcval[1] < 0))
        {
            tha = thc;
            thc = (tha + thb) / 2;
        }

    } // End the bisections.


    // Compute the final stats.
    InfoVal = 0;

    // Loop over items.
    for (j = 0; j < J; j++)
    {
        // Step probabilities.
        a = ipar(j, 0);
        P[0] = 1;
        Q[0] = 0;

        for (m = 0; m < M; m++)
        {
            b = ipar(j, (m + 1));
            logit = a * (thc - b);
            P[m + 1] = exp(logit) / (1 + exp(logit));
            Q[m + 1] = 1 - P[m + 1];
        }

        P[M + 1] = 0;
        Q[M + 1] = 1;

        for (m = 0; m < (M + 1); m++)
        {
            Pr = P[m] - P[(m + 1)];
            D1 = a * (P[m] * (1 - P[m]) - P[m + 1] * (1 - P[m + 1]));
            InfoVal = InfoVal + pow(D1, 2) / Pr;
        }
    } // End loop over items.


    // Prepare the return.
    thHat[0] = thc;
    thHat[1] = 1 / sqrt(InfoVal);

    return (thHat);
}

// endregion



// region Theta WML and PolyLzStar (wrappers)

Rcpp::NumericMatrix Statistics::FullThetaWml(const Rcpp::NumericMatrix &data, const Rcpp::NumericMatrix &item_parameters)
{
    int rows = data.nrow();
    Rcpp::NumericMatrix theta_wml_stats(rows, 4);


    for (int row = 0; row < rows; ++row)
    {
        // Estimate the theta by WML for a single response vector.
        Rcpp::NumericVector estimated_theta = EstimateThetaWml(data(row, Rcpp::_), item_parameters);


        // Perform some checks on the estimated theta vector to determine
        // whether the bisection method must be applied.
        if (estimated_theta(2) == 0 || estimated_theta(3) == 0)
        {
            // If execution hits this point it means that the estimation by means
            // of Newton-Raphson optimization did not converge. Thus, we fall
            // back on the bisection method and we flag the last two
            // columns with value 2, so we make this visible.
            Rcpp::NumericVector bisection = Statistics::EstimateThetaBisection(data(row, Rcpp::_), item_parameters);
            theta_wml_stats(row, Rcpp::_) = Rcpp::NumericVector::create(bisection(0), bisection(1), 2, 2);
        }
        else
        {
            theta_wml_stats(row, Rcpp::_) = estimated_theta;
        }
    }

    return theta_wml_stats;
}



Rcpp::NumericMatrix Statistics::FullPolyLzStar(const Rcpp::NumericMatrix &data, const Rcpp::NumericMatrix &item_parameters, const Rcpp::DoubleVector &theta)
{
    int rows = data.nrow();

    Rcpp::NumericMatrix poly_lz_star(rows, 2);

    for (int row = 0; row < rows; ++row)
    {
        poly_lz_star(row, Rcpp::_) = EstimatePolyLzStar(data(row, Rcpp::_), item_parameters, theta(row));
    }

    return poly_lz_star;
}

// endregion



// Cohen's d effect size statistic.
double Statistics::CohenEffectSize(const Rcpp::NumericVector &firstVector, const Rcpp::NumericVector &secondVector)
{
    return (Rcpp::mean(firstVector) - Rcpp::mean(secondVector)) / sqrt((pow(Rcpp::sd(firstVector), 2) + pow(Rcpp::sd(secondVector), 2)) / 2);
}