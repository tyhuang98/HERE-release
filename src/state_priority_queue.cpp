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

#include "state_priority_queue.h"


namespace HERE
{

    template <class SSR, typename Scalar>
    StatePriorityQueue<SSR, Scalar>::StatePriorityQueue(OptProblem op):
            optProblem(op),
            head(NULL), tail(NULL), m_size(0) {}

    template <class SSR, typename Scalar>
    StatePriorityQueue<SSR, Scalar>::~StatePriorityQueue()
    {
        Node *node, *tmp;

        node = head;
        while(node)
        {
            tmp = node;
            node = node->right;
            delete tmp;
        }
    }


    template <class SSR, typename Scalar>
    HERE::utils::SearchState<SSR, Scalar> *StatePriorityQueue<SSR, Scalar>::pop()
    {
        HERE::utils::SearchState<SSR, Scalar> *res;

        if( head==NULL )
        {
            return NULL;
        }

        Node *tmp=head;
        res = head->state;

        head = head->right;
        if (head != NULL)
        {
            head->left = NULL;
        }
        else
        {
            tail = NULL;
        }

        tmp->state=NULL;
        delete tmp;
        m_size--;
        return res;
    }


    template <class SSR, typename Scalar>
    void StatePriorityQueue<SSR, Scalar>::push(HERE::utils::SearchState<SSR, Scalar> *state)
    {
        Node *node, *new_node;
        new_node = new Node(state);

        if (m_size==0)
        {
            head = tail = new_node;
            m_size = 1;
            return;
        }

        // Traverse queue from tail to insert node in descent order
        node = tail;
        while (node)
        {
            if ((optProblem==MAXIMISATION && node->state->bnd >= state->bnd) ||
                (optProblem==MINIMISATION && node->state->bnd <= state->bnd) )
            {
                // Inserting at the end of the queue
                if (node==tail)
                {
                    tail = new_node;
                }
                    // Inserting node between head and tail nodes
                else
                {
                    new_node->right = node->right;
                    new_node->right->left = new_node;
                }

                // Link new_node to the left node
                new_node->left = node;
                node->right = new_node;

                m_size++;
                return;
            }
            node = node->left;
        }

        // Inserting at the head
        new_node->right = head;
        new_node->right->left = new_node;
        head = new_node;
        m_size++;
    }

    template <class SSR, typename Scalar>
    void StatePriorityQueue<SSR, Scalar>::prune(int curbest)
    {
        Node *node, *tmp_node;

        // Traverse queue from tail to delete nodes
        node = tail;
        while (node)
        {
            // Stop condition
            if ( (optProblem==MAXIMISATION && node->state->bnd > curbest) ||
                 (optProblem==MINIMISATION && node->state->bnd < curbest))
            {
                tail = node;
                tail->right = NULL;
                return;
            }

            // Delete visited node
            m_size--;
            tmp_node = node;
            node = node->left;
            delete tmp_node;
        }

        // Boundary case: All nodes were deleted
        head = tail = NULL;
        assert(m_size>=0);
    }

    template <class SSR, typename Scalar>
    unsigned int StatePriorityQueue<SSR, Scalar>::size() const
    {
        assert(m_size>=0);
        return m_size;
    }

} // End namespace reg

