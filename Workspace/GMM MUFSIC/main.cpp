#include <armadillo>
#include <iostream>

int main (void) {
	std::cout << "Hello Project" << std::endl;
	arma::mat a = arma::randu<arma::mat>(5, 5);
	arma::vec b = arma::randu<arma::vec>(5);
	arma::vec x;
	arma::solve(x, a, b);
	std::cout << x << std::endl;
	return 0;
}
