#include <iostream>
#include <boost/numeric/interval.hpp>

int main()
{
  boost::numeric::interval<double> i;
  std::cout << "Hello world " << i.lower() << ", " << i.upper();
}
