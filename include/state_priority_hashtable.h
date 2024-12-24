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

#ifndef REG_STATE_PRIORITY_HASHTABLE_H_
#define REG_STATE_PRIORITY_HASHTABLE_H_

#include "utils.h"
#include <cstddef>

#ifdef __linux__ //asuming gcc
#include <tr1/unordered_map>
using namespace std::tr1;
#else
#include <unordered_map> //asuming clang
#endif

namespace HERE {

    template <class SSR, typename Scalar=int >
    class Queue
    {
    private:
        class Node
        {
        public:
            Node *next;
            HERE::utils::SearchState<SSR, Scalar> *ss;

            Node(HERE::utils::SearchState<SSR, Scalar> *ss): ss(ss){}
            ~Node()
            {
                delete ss;
            }

        };
    public:
        Queue();
        ~Queue();

        Node *head;
        Node *tail;

        int size;
        HERE::utils::SearchState<SSR, Scalar> *pop();
        void push(HERE::utils::SearchState<SSR, Scalar> *ss);
        void stack(HERE::utils::SearchState<SSR, Scalar> *ss); //not the best design...

        void dump(std::ofstream& ofs) const;
    }; // Queue


    struct Hash
    {
        std::size_t operator()(const int& q) const
        {
            unsigned int key=q;
            key = ((key >> 16) ^ key) * 0x45d9f3b;
            key = ((key >> 16) ^ key) * 0x45d9f3b;
            key = ((key >> 16) ^ key);
            return key;
        }
    };

//---------------------------------------------
//      State Priority Hashtable
//---------------------------------------------
    template <class SSR, typename Scalar=int >
    class StatePriorityHashtable
    {

    public:
        explicit StatePriorityHashtable(int bucketSize);
        ~StatePriorityHashtable();

        HERE::utils::SearchState<SSR, Scalar> *pop();
        void push(HERE::utils::SearchState<SSR, Scalar> *state);
        void prune(int curbest);
        unsigned int size() const;

        //dump table SSR to file
        void dump(std::ofstream& ofs) const;

    private:

#ifdef __linux__ //asuming gcc
        typedef std::tr1::unordered_map<Scalar, Queue<SSR, Scalar>*, Hash> Map;
#else
        typedef std::unordered_map<Scalar, Queue*, Hash> Map;
#endif

        int m_max_upbnd;
        Queue<SSR, Scalar> *m_max_upbnd_queue;
        unsigned int m_size; // number of states

        Map m_map;
    };

} // End namespace reg

#endif
