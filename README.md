# Monte Carlo simulation for stock price dynamics using a recursive form of the Geometric Brownian Motion SDE
The aim of this project is to develop a quantitative model that can measure and visualize the inherent uncertainty and randomness in the movement of stock prices. By analyzing the fluctuations in the market, the model provides valuable insights that can assist in evaluating the behavior of securities and forming well-supported hypotheses around future actions that could be taken.

It should be emphasized that the model is not intended to predict the future movement of stock prices with complete accuracy, but rather to provide probabilistic estimates that can inform decision-making. As with any quantitative model, the accuracy of the output is subject to the quality and quantity of the input data, as well as the assumptions and limitations of the model itself. Therefore, it is important to approach the model's outputs with a healthy degree of skepticism and critical thinking, and to use it in combination with other sources of information and expertise.

## References:
* Geometric brownian motion. (n.d.). https://www.quantstart.com/articles/Geometric-Brownian-Motion/
* Lewinson, E. (2022). Python for finance cookbook: Over 80 powerful recipes for effective financial data analysis. Packt Publishing.
* K. (2013). Simulating brownian motion (bm) and Geometric Brownian motion (GBM). http://www.columbia.edu/~ks20/4404-Sigman/4404-Notes-sim-BM.pdf
* Rouah, F. (2010). Euler and Milstein discretization - frouah.com. http://www.frouah.com/finance%20notes/Euler%20and%20Milstein%20Discretization.pdf
* MIT OpenCourseWare. (2013). Lecture 21: Stochastic differential equations: Topics in mathematics with applications in finance: Mathematics. MIT OpenCourseWare. https://ocw.mit.edu/courses/18-s096-topics-in-mathematics-with-applications-in-finance-fall-2013/resources/lecture-21-stochastic-differential-equations/ 
* Shreve, S. E. (2004). Stochastic Calculus for Finance II: Continuous-Time Models. Springer.
* Wilmott, P., Howison, S., & Dewynne, J. N. (2013). The Mathematics of Financial Derivatives: A Student Introduction. Cambridge University Press.
* Hull, J. C. (2017). Options, Futures, and Other Derivatives (10th ed.). Pearson Education.

## A breif description of Geometric brownian motion and the derived recursive form used in this model for estimating geometric borwnian motion in stock price path dynamics:
### Geometric brownian motion:

Geometric Brownian Motion is a continuous time stochastic processs used to describe the stochastic movement of stock prices. 

Stock prices are considered to be a stochastic process beause they are subject to random fluctuations that are influenced by large number of uncontrollable factors (i.e. economic new, company performance, and investor sentiment). 

As with any quantitive model used for informing investment decisions there are limitations. Of note, Geometric Brownian Motion Model has to following limiting assumptions:
* Volatility of the price process is constant over time. 
* Lack of mean reversion. The GBM model assumes that asset prices follow a random walk, which means that they move randomly in the absence of any other information.
* Normality assumption. Which is in conflict with decades of empirical evidence that suggests asset returns often exhibit non-normal distributions, such as fat tails or skewness.
* Lack of jumps. Assuming that prices change continuously over time with no sudden jumps or shocks caused by stock market crashes, flash crashes, or news-driven price spikes

That being said even after considering the limitations of the model, it continues to provide a well defined starting point for modeling asset prices.

### The GBM model is expressed by the following stochastic differential equation:
<p align="center">
$$dS(t) = μS(t)dt + σS(t)dW(t)$$
</p>

where:
* $S(t)$ is the price of the asset at time t
* $μ$, (drift), is the expected rate of return on the asset
* $σ$, (diffusion), is the volatility of the asset
* $W(t)$ is a Wiener process, which is a mathematical representation of Brownian motion

In this case W(t) is simulated using using a simple Monte Carlo simulation by generating a series normally distrubted random numbers which are used to update the state of the stochastic process at each time step. Although not an exact weiner process, the simulation will have similar properties such as a continuous and unpredicatable path at each time step.

### The derivation of the recursive form of the GBM SDE (shown below) using Euler-Maruyama's Discretization:
<p align="center">
$$S(t+dt)=S(t)\cdot\exp\left((\mu-\frac{1}{2}\sigma^2)dt+\sigma\sqrt{dt}\cdot Z\right)$$
</p>

where:
* $S(t)$ is the stock price at time $t$
* $dt$ is the time step size
* $μ$ is the drift or expected stock return
* $σ$ is the volatility or diffusion
* $Z$ is a standard normal random variable
* $exp$ is the exponential function

and:
* $(\mu-\frac{1}{2}\sigma^2)dt$ represents the expected change in the log price over the time step

and:
* $\sigma\sqrt{dt}\cdot Z$  represents the random component of the log price change.

### Deriving the Recursive Form of Geometric Brownian Motion using the Euler-Maruyama Discretization method:

We start with the original stochastic differential equation $$dS(t) = μS(t)dt + σS(t)dW(t)$$
where (as mentioned above) $S(t)$ is the stock price at time $t$, $\mu$ is the drift rate (the expected rate of return), $\sigma$ is the volatility (the standard deviation of the rate of return), $W(t)$ is a standard Brownian motion (a stochastic process with independent and normally distributed increments), and $dt$ is an infinitesimal time step.

To obtain the recursive form, we need to solve the SDE. One way to do this is to use the Euler-Maruyama method, which involves discretizing time and approximating the stochastic differential equation by a difference equation. Note: A difference equation is a mathematical equation that describes a sequence of values in terms of the values that precede it.

Below are the steps of Euler-Maruyama Discretization:

1. Discretize time into small intervals of length dt. We want to find the value of $S$ at time $t+dt$, given the value of $S$ at time $t$. We can write: $$S(t+dt) = S(t) + dS(t)$$
2. Substitute the SDE for dS(t): $$S(t+dt) = S(t) + μS(t)dt + σS(t)dW(t)$$
3. Approximate the increment $dW(t)$ using a normal distribution with mean $0$ and variance $dt$: $$dW(t) \approx N(0,dt)$$
4. Define a standard normal random variable $Z$ as:
$$Z \approx N(0,1)$$
5. Rewrite the equation using $dW(t) = sqrt(dt)\cdot Z$: $$S(t+dt) = S(t) + μS(t)dt + σS(t)\sqrt{dt}\cdot Z$$
6. Use the fact that $Z$ is normally distributed to obtain: $$S(t+dt) = S(t) + μS(t)dt + σS(t)\sqrt{dt}\cdot Z \approx S(t)\cdot\left(1 + μdt + σ\sqrt{dt}\cdot Z\right)$$
7. Finally, exponentiate both sides of the equation to obtain the recursive form of the GBM: $$S(t+dt)=S(t)\cdot\exp\left((\mu-\frac{1}{2}\sigma^2)dt+\sigma\sqrt{dt}\cdot Z\right)$$

The final equation is the discrete approximation to the original GBM SDE, using the Euler-Maruyama method. It tells us that the predicted stock price after a small time interval $dt$, $S(t + dt)$, can be obtained by multiplying the current stock price $S(t)$ by a random factor $\exp \left( (\mu - \frac{1}{2}\sigma^2) dt + \sigma \sqrt{dt} Z \right)$, where $Z$ is a standard normal random variable. This random factor captures the stochastic fluctuations in the stock price due to the Wiener process $dW(t)$.

## Setup:
* Import necessary packages (see requirements file).
* Define the security for which you would like to simulate price dynamics for. 
* Choose the time frame from which you will train and test the Monte Carlo Simulation model.
* Download you data using the yf.download() function.
* specify the Adj Close column to account for an income accrual events in the price of the security (dividends, short term and long term capital gains, corp actions, etc.). adj_close = x['Adj Close']
* Calculate the daily percentage change in adjusted close prices and drop the first row which contains null values. returns = adj_close.pct_change().dropna()

### Preparing your data:
* Split your returns dataframe into your training set and test set
* Set up your Monte Carlo simulation parametes: Including the length of the test period (T), the number of time steps (N), the initial stock price (S_0), the number of simulations (N_SIM), and the mean (mu) and standard deviation (sigma) of the returns during the training period.

### The simulate_gbm function:
* This function simulates the Geometric Brownian Motion (GBM) of the selcted security using a Monte Carlo simulation. The function input parameters are :
* s_0: the initial stock price
* mu: the expected return of the stock
* sigma: the volatility of the stock
* n_sims: the number of simulations to run
* T: the time horizon
* N: the number of time steps to divide the horizon into

The function calculates the time step dt by dividing T by N. It then generates an array of random normal variables dW with a scale of the square root of dt using the NumPy np.random.normal function.

Next, the simulate_gbm function uses NumPy's np.cumsum function to compute the cumulative sum of dW along the second axis to obtain an array W of Brownian motion values. It also creates a matrix time_steps consisting of repeated copies of the time steps.

Finally, the function computes the simulated stock prices S_t at each time step by applying the GBM formula using the np.exp and np.insert functions to insert the initial stock price s_0 at the beginning of each row of the array. The function returns the simulated stock prices S_t.

### Refining and Vizualizing your results against to original mean return:
* LAST_TRAIN_DATE = train.index[-1].date() - sets the variable LAST_TRAIN_DATE to the date of the last entry in the train DataFrame's index. train is a DataFrame containing historical data up until the last date of the training set.
* FIRST_TEST_DATE = test.index[0].date() - sets the variable FIRST_TEST_DATE to the date of the first entry in the test DataFrame's index. test is a DataFrame containing historical data starting from the first date of the testing set.
* LAST_TEST_DATE = test.index[-1].date() - sets the variable LAST_TEST_DATE to the date of the last entry in the test DataFrame's index.
* PLOT_TITLE = (f'{ticker} Simulation 'f'({FIRST_TEST_DATE}:{LAST_TEST_DATE})') - creates a string that will be used as the title of the plot. 
* selected_indices selects the date range for which we have both actual stock price data and simulated data.
* gbm_simulations_df creates a pandas DataFrame of the GBM simulation results, with the date range from selected_indices as the index.
* The plot function is then used to plot the DataFrame gbm_simulations_df with a specified figsize and an alpha of 0.2 (to make the simulation lines more transparent), and with no legend (legend=False).
* line_1 and line_2 are created using the plot function, with line_1 representing the simulated mean of the GBM simulations, and line_2 representing the actual stock price path.
* The set_title function sets the plot title to PLOT_TITLE.
* The legend function adds a legend to the plot, with line_1 labeled as 'Simulated Mean' and line_2 labeled as 'Actual Price Path'.
* The plot is then displayed using plt.show().
