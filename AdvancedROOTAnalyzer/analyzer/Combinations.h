/*! \file Header file containing the combination info necessary for alpha_T calculations.
 */

#ifndef Combinations_h
#define Combinations_h

#include <sstream>
#include <iostream>
#include <vector>

namespace Combinations {

    //UInt_t nmax = 10; // Determined by how much is coded! A move up to 11 would require double digits... or chars...
    
    inline int factorial(int n) {
        int fact = 1;
        for (int i=n; i>=1; i--) fact = fact * i;
        return fact;
    }

    // Non recursive template function
    // Algorithm taken from - thanks! :-)
    // http://www.codeguru.com/cpp/cpp/algorithms/combinations/article.php/c5117/
    // http://www.codeproject.com/KB/recipes/combinations_cpp_2.aspx
    template <class BidIt>

    inline bool next_combination(BidIt n_begin, BidIt n_end, BidIt r_begin, BidIt r_end) {
        
        bool boolmarked=false; // ??? What does this mean?
        
        BidIt r_marked;
        
        BidIt n_it1 = n_end; // The end of the n collection
        --n_it1; // Go one back?
        
        BidIt tmp_r_end=r_end; // The end of the r collection
        --tmp_r_end;
  
        for(BidIt r_it1=tmp_r_end; r_it1!=r_begin || r_it1==r_begin; --r_it1,--n_it1) {
            if(*r_it1==*n_it1 ) {
	            if(r_it1!=r_begin) {//to ensure not at the start of r sequence
		            boolmarked=true;
		            r_marked=(--r_it1);
		            ++r_it1;//add it back again
		            continue;
		        }
	            else { // it means it is at the start the sequence, so return false
		            return false;
                }
	        }
	        else //if(*r_it1!=*n_it1 )
	        {
	            //marked code
	            if(boolmarked==true)
		        {
		            //for loop to find which marked is in the first sequence
		            BidIt n_marked;//mark in first sequence
		            for (BidIt n_it2=n_begin;n_it2!=n_end;++n_it2) {
		                if(*r_marked==*n_it2) {n_marked=n_it2;break;}
                    }
		            BidIt n_it3=++n_marked;
		            for  (BidIt r_it2=r_marked;r_it2!=r_end;++r_it2,++n_it3) {
		                *r_it2=*n_it3;
		            }
		            return true;
		        }
	            for(BidIt n_it4=n_begin; n_it4!=n_end; ++n_it4) {
		            if(*r_it1==*n_it4) {
		                *r_it1=*(++n_it4);
		                return true;
		            }
	            }
            }
	    }

        return true;//will never reach here
    } // end of next_combination method

    inline void mycombinations( UInt_t size, 
                                std::vector< std::vector<UInt_t> >& combo1, 
		                        std::vector< std::vector<UInt_t> >& combo2  ) {
        size = size < 10 ? size : 10;
  
        combo1.clear();
        combo2.clear();
  
        if ( size == 0 ) { return; }
        if ( size == 1 ) { 
            combo1.push_back( std::vector<UInt_t>(1,0) );
            combo2.push_back( std::vector<UInt_t>() );
            return; 
        }
  
        // Fill the list of indices for the n objects and copy to array
        std::stringstream ss1;
        for ( UInt_t ii = 0; ii < size; ++ii ) { ss1 << ii; }
        char cc1[100];
        strcpy( cc1, ss1.str().c_str() );
    
        char cc2[100];

        for ( UInt_t j = 0; j < size/2; ++j) {
            std::stringstream ss2;
            for ( UInt_t jj = 0; jj < j+1; ++jj) {
                //for ( UInt_t jj = 0; jj < size/2; ++jj ) { ss2 << jj; }
                ss2 << jj;
            }
            strcpy( cc2, ss2.str().c_str() );

            UInt_t size1 = ss1.str().size();
            UInt_t size2 = ss2.str().size();
  
            do { 
                std::vector<UInt_t> tmp1;
                for ( UInt_t ii = 0; ii < size2; ++ii ) { 
                    std::stringstream sss; sss << cc2[ii];
                    UInt_t num; sss >> num;
                    tmp1.push_back( num ); 
                }
                combo1.push_back( tmp1 );
                std::vector<UInt_t> tmp2;
                for ( UInt_t ii = 0; ii < size1; ++ii ) { 
                    if ( std::find( tmp1.begin(), tmp1.end(), ii ) == tmp1.end() ) { 
	                    tmp2.push_back(ii); 
                    }
                }
                combo2.push_back( tmp2 );
            } while ( next_combination( cc1, cc1 + size1, cc2, cc2 + size2 ) );
  
        } //end of loop over r values
 
        //if ((size>0)&&(float(size)/2.0-float(int(size/2))==0.)) {
        if ((size>0)&&((size % 2)==0)) {
            UInt_t numpop = factorial(size) / (factorial(size/2)*factorial(size/2));
            //std::cout << numpop/2 << std::endl;
            for (UInt_t i=0;i<numpop/2;i++) {
                combo1.pop_back();
                combo2.pop_back();
            }
        }
  
        if (0) //( edm::isDebugEnabled() ) 
        {
            std::stringstream ss;
            ss << " size=" << size << std::endl;
            UInt_t combo_size = combo1.size() < combo2.size() ? combo1.size() : combo2.size();
            ss << " combo1.size()=" << combo1.size() 
               << " combo2.size()=" << combo2.size()
               << std::endl;
            for ( UInt_t ii = 0; ii < combo_size; ++ii ) { 
                ss << " combo1[" << ii << "].size()=" << combo1[ii].size()
	               << ", values=";
                for ( UInt_t jj = 0; jj < combo1[ii].size(); ++jj ) { 
	                ss << combo1[ii][jj] << ",";
                }
                ss << " combo2[" << ii << "].size()=" << combo2[ii].size()
	               << ", values=";
                for ( UInt_t jj = 0; jj < combo2[ii].size(); ++jj ) { 
	                ss << combo2[ii][jj] << ",";
                }
                ss << std::endl;
            }
            std::cout << ss.str() << std::endl;
            //LogTrace("TEST") << ss.str();
        } // end of debug enabled check
  
    } // end of mycombinations



/* Old, hard-coded method

  // n =           0  1  2  3  4   5   6
  UInt_t rmax[] = {0, 0, 1, 1, 2,  2,  3};
  UInt_t jmax[] = {0, 0, 1, 3, 4, 10, 15};

  int la2[][1][1] = {{{0}}};
  int lb2[][1][1] = {{{1}}};

  int la3[][3][2] = {{{0, -1}, {1, -1}, {2, -1}}};
  int lb3[][3][2] = {{{1,  2}, {0,  2}, {0,  1}}};

  int la4[][4][3] = {{{ 0,-1,-1}, { 1,-1,-1}, { 2,-1,-1}, { 3,-1,-1}},
		     {{ 0, 1,-1}, { 0, 2,-1}, { 0, 3,-1}, {-1,-1,-1}}};
  int lb4[][4][3] = {{{ 1, 2, 3}, { 0, 2, 3}, { 0, 1, 3}, { 0, 1, 2}},
		     {{ 2, 3,-1}, { 1, 3,-1}, { 1, 2,-1}, {-1,-1,-1}}};

  int la5[][10][4] = {{{ 0,-1,-1,-1}, { 1,-1,-1,-1}, { 2,-1,-1,-1},
		       { 3,-1,-1,-1}, { 4,-1,-1,-1},
		       {-1,-1,-1,-1}, {-1,-1,-1,-1}, {-1,-1,-1,-1},
		       {-1,-1,-1,-1}, {-1,-1,-1,-1}},
		      {{ 0, 1,-1,-1}, { 0, 2,-1,-1}, { 0, 3,-1,-1},
		       { 0, 4,-1,-1}, { 1, 2,-1,-1},
		       { 1, 3,-1,-1}, { 1, 4,-1,-1}, { 2, 3,-1,-1},
		       { 2, 4,-1,-1}, { 3, 4,-1,-1}}};
  
  int lb5[][10][4] = {{{ 1, 2, 3, 4}, { 0, 2, 3, 4}, { 0, 1, 3, 4},
		       { 0, 1, 2, 4}, { 0, 1, 2, 3},
		       {-1,-1,-1,-1}, {-1,-1,-1,-1}, {-1,-1,-1,-1},
		       {-1,-1,-1,-1}, {-1,-1,-1,-1}},
		      {{ 2, 3, 4,-1}, { 1, 3, 4,-1}, { 1, 2, 4,-1},
		       { 1, 2, 3,-1}, { 0, 3, 4,-1},
		       { 0, 2, 4,-1}, { 0, 2, 3,-1}, { 0, 1, 4,-1},
		       { 0, 1, 3,-1}, { 0, 1, 2,-1}}};

  int la6[][15][5] = {{{ 0,-1,-1,-1,-1},{ 1,-1,-1,-1,-1},{ 2,-1,-1,-1,-1},
		       { 3,-1,-1,-1,-1},{ 4,-1,-1,-1,-1},{ 5,-1,-1,-1,-1},
		       {-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},
		       {-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},
		       {-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1}},
		      {{ 0, 1,-1,-1,-1},{ 0, 2,-1,-1,-1},{ 0, 3,-1,-1,-1},
		       { 0, 4,-1,-1,-1},{ 0, 5,-1,-1,-1},{ 1, 2,-1,-1,-1},
		       { 1, 3,-1,-1,-1},{ 1, 4,-1,-1,-1},{ 1, 5,-1,-1,-1},
		       { 2, 3,-1,-1,-1},{ 2, 4,-1,-1,-1},{ 2, 5,-1,-1,-1},
		       { 3, 4,-1,-1,-1},{ 3, 5,-1,-1,-1},{ 4, 5,-1,-1,-1}},
		      {{ 0, 1, 2,-1,-1},{ 0, 1, 3,-1,-1},{ 0, 1, 4,-1,-1},
		       { 0, 1, 5,-1,-1},{ 0, 2, 3,-1,-1},{ 0, 2, 4,-1,-1},
		       { 0, 2, 5,-1,-1},{ 0, 3, 4,-1,-1},{ 0, 3, 5,-1,-1},
		       { 0, 4, 5,-1,-1},{-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},
		       {-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1}}};
  int lb6[][15][5] = {{{ 1, 2, 3, 4, 5},{ 0, 2, 3, 4, 5},{ 0, 1, 3, 4, 5},
		       { 0, 1, 2, 4, 5},{ 0, 1, 2, 3, 5},{ 0, 1, 2, 3, 4},
		       {-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},
		       {-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},
		       {-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1}},
		      {{ 2, 3, 4, 5,-1},{ 1, 3, 4, 5,-1},{ 1, 2, 4, 5,-1},
		       { 1, 2, 3, 5,-1},{ 1, 2, 3, 4,-1},{ 0, 3, 4, 5,-1},
		       { 0, 2, 4, 5,-1},{ 0, 2, 3, 5,-1},{ 0, 2, 3, 4,-1},
		       { 0, 1, 4, 5,-1},{ 0, 1, 3, 5,-1},{ 0, 1, 3, 4,-1},
		       { 0, 1, 2, 5,-1},{ 0, 1, 2, 4,-1},{ 0, 1, 2, 3,-1}},
		      {{ 3, 4, 5,-1,-1},{ 2, 4, 5,-1,-1},{ 2, 3, 5,-1,-1},
		       { 2, 3, 4,-1,-1},{ 1, 4, 5,-1,-1},{ 1, 3, 5,-1,-1},
		       { 1, 3, 4,-1,-1},{ 1, 2, 5,-1,-1},{ 1, 2, 4,-1,-1},
		       { 1, 2, 3,-1,-1},{-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},
		       {-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1},{-1,-1,-1,-1,-1}}};

*/

} // ~namespace Combinations


#endif // ~Combinations_h
