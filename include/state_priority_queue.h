/**
*                    GUARANTEED OUTLIER REMOVAL
*         FOR POINT CLOUD REGISTRATION WITH CORRESPONDENCES
*
*
*
* Copyright (C) 2016 Alvaro PARRA BUSTOS (aparra at cs.adelaide.edu.au)
* School of Computer Science, The University of Adelaide, Australia
* The Australian Center for Visual Technologies
* http://cs.adelaide.edu.au/~aparra
*
* The source code, binaries and demo is distributed for academic use only,
* for any other use contact the authors.
*/

#ifndef REG_STATE_PRIORITY_QUEUE_H_
#define REG_STATE_PRIORITY_QUEUE_H_

#include "utils.h"
#include <cstddef> //NULL

#include <time.h>
#include <stdlib.h>
#include <cstdio>
#include <assert.h>

namespace HERE {

    template <class SSR, typename Scalar=int >
    class StatePriorityQueue
    {
    public:
        enum OptProblem{MINIMISATION, MAXIMISATION};

    private:

        class Node
        {
        public:
            HERE::utils::SearchState<SSR, Scalar> *state;
            Node *left, *right;

            Node(): state(NULL), left(NULL), right(NULL) {}
            Node(HERE::utils::SearchState<SSR, Scalar> *state): state(state), left(NULL), right(NULL){}
            ~Node() {if (state!=NULL) delete state;}
        };

        const OptProblem optProblem;
        Node *head, *tail;
        unsigned int m_size;

    public:
        StatePriorityQueue(OptProblem op=MAXIMISATION);
        ~StatePriorityQueue();

        HERE::utils::SearchState<SSR, Scalar> *pop();
        void push(HERE::utils::SearchState<SSR, Scalar> *state);

        /**
         * @brief Remove and free states with upper bound lower or equal to lwbnd.
         * @param lwbnd Known lower bound.
         */
        void prune(int curbest);

        unsigned int size() const;
    };

} // End namespace reg

#endif

