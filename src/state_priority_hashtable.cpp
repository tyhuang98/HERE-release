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
#include "state_priority_hashtable.h"

namespace HERE {

    template <class SSR, typename Scalar >
    Queue<SSR, Scalar>::Queue(): size(0) {}

    template <class SSR, typename Scalar >
    Queue<SSR, Scalar>::~Queue() {
        Node *node, *tmp;
        node = head;
        for (int i = 0; i < size; i++) {
            tmp = node;
            node = node->next;
            delete tmp;
        }
    }

    template<class SSR, typename Scalar>
    utils::SearchState<SSR, Scalar> *Queue<SSR, Scalar>::pop() {

//    mxAssert(size>0, "pop in an empty queue");
        assert(size > 0);

        Node *tmp = head;
        HERE::utils::SearchState<SSR, Scalar> *ss = head->ss;
        head = head->next;
        tmp->ss = NULL;
        delete tmp;
        size--;

//    mxAssert(size>=0, "Invalid size");
        assert(size >= 0);
        return ss;

    }

    template<class SSR, typename Scalar>
    void Queue<SSR, Scalar>::push(utils::SearchState<SSR, Scalar> *ss) {
        Node *node;

        if (size == 0) {
            head = tail = new Node(ss);
            size = 1;
            return;
        }

        node = tail;
        // append at the end
        tail = new Node(ss);
        node->next = tail;
        size++;
    }

    template<class SSR, typename Scalar>
    void Queue<SSR, Scalar>::stack(utils::SearchState<SSR, Scalar> *ss) {
        Node *node;

        if (size == 0) {
            head = tail = new Node(ss);
            size = 1;
            return;
        }

        node = head;
        // stack at the begining
        head = new Node(ss);
        //node->next = tail;
        head->next = node;
        size++;
    }

    template<class SSR, typename Scalar>
    void Queue<SSR, Scalar>::dump(std::ofstream &ofs) const {
        Node *node = head;
        HERE::utils::SearchState<SSR, Scalar> *ss;

        for (int i = 0; i < size; i++) {
            ss = node->ss;
//            ofs << *ss << "\n";
            node = node->next;
        }
    }



//-------------------------------
//     State Priority Hashtable
//-------------------------------

    template <class SSR, typename Scalar >
    StatePriorityHashtable<SSR, Scalar>::StatePriorityHashtable(int bucketSize):
            m_max_upbnd(0), m_size(0) {
#ifdef __linux__ //asuming gcc
        m_map = std::tr1::unordered_map<Scalar, Queue<SSR, Scalar> *, Hash>(bucketSize, Hash());
#else
        m_map = std::unordered_map<Scalar, Queue*, Hash>(bucketSize, Hash());
#endif
    }

    template <class SSR, typename Scalar >
    StatePriorityHashtable<SSR, Scalar>::~StatePriorityHashtable() {
        typename Map::iterator it;
        for (it = m_map.begin(); it != m_map.end(); ++it) {
            delete it->second;
        }
    }


    template <class SSR, typename Scalar >
    HERE::utils::SearchState<SSR, Scalar> *StatePriorityHashtable<SSR, Scalar>::pop() {
        HERE::utils::SearchState<SSR, Scalar> *searchState;
        Queue<SSR, Scalar> *queue;
        //SSR ssr;

//    mxAssert(m_size>0 , "pop in an Empty table");
//    mxAssert(m_map.count(m_max_upbnd)>0, "No queue in table");
//    mxAssert(m_map[m_max_upbnd]==m_max_upbnd_queue, "m_max_upbnd_queue");
        assert(m_size > 0);
        assert(m_map.count(m_max_upbnd) > 0);
        assert(m_map[m_max_upbnd] == m_max_upbnd_queue);

        queue = m_max_upbnd_queue;
        searchState = queue->pop();
        m_size--;

        // Remove queue if it is empty
        if (queue->size == 0) {
            m_map.erase(m_max_upbnd);
            delete queue;

            // Update max_upbnd
            if (m_size) {
                m_max_upbnd--;
                while (m_map.count(m_max_upbnd) == 0) {
                    m_max_upbnd--;
                }

                m_max_upbnd_queue = m_map[m_max_upbnd];
            } else {
                m_max_upbnd = 0;
            }
        }

        return searchState;
    }

    template <class SSR, typename Scalar >
    void StatePriorityHashtable<SSR, Scalar>::push(HERE::utils::SearchState<SSR, Scalar> *state) {
//        mxAssert(state->bnd >= 0, "inserting state with upbnd<=0");
//        static_assert(state->bnd >= 0, "inserting state with upbnd<=0");
        Queue<SSR, Scalar> *q;

        if (m_map.count(state->bnd) == 0) {
            q = new Queue<SSR, Scalar>();
            m_map[state->bnd] = q;
        } else {
            q = m_map[state->bnd];
//                mxAssert(q->size > 0, "Invalid queue size in hashtable");
//                mxAssert(m_map.count(state->bnd) > 0, "Map error");
//            static_assert(q->size > 0, "Invalid queue size in hashtable");
//            static_assert(m_map.count(state->bnd) > 0, "Map error");
        }

        q->push(state);
//            mxAssert(q->size > 0, "Invalid queue size after push");
//        static_assert(q->size > 0, "Invalid queue size after push");

        // Update max_upbnd
        if (state->bnd > m_max_upbnd) {
            m_max_upbnd = state->bnd;
            m_max_upbnd_queue = q;
        }

        m_size++;
    }

    template <class SSR, typename Scalar >
    void StatePriorityHashtable<SSR, Scalar>::prune(int lwbnd) {
        Queue<SSR, Scalar> *queue;
        typename Map::iterator it;

        it = m_map.begin();

        while (it != m_map.end()) {
            if (it->first <= lwbnd) {
                queue = it->second;
                // remove and delete queue
                m_size -= queue->size;
                delete queue;

                it = m_map.erase(it);
            } else {
                it++;
            }
        }

        // Update max_upbnd
        if (m_size == 0) {
            m_max_upbnd = 0;
        }
    }


    template <class SSR, typename Scalar >
    void StatePriorityHashtable<SSR, Scalar>::dump(std::ofstream &ofs) const {
        Queue<SSR, Scalar> *queue;
        typename Map::const_iterator it;

        it = m_map.begin();

        while (it != m_map.end()) {
            queue = it->second;
            //queue->dump();
            queue->dump(ofs);
            it++;
        }
    }

    template <class SSR, typename Scalar >
    unsigned int StatePriorityHashtable<SSR, Scalar>::size() const {
        return m_size;
    }

    template class StatePriorityHashtable<utils::state1Dof, int>;

}
