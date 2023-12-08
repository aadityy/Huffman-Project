
#include "pq.h"

#include <stdbool.h>
#include <stdio.h>
/* put in pq.c */
typedef struct ListElement ListElement;
struct ListElement {
    Node *tree;
    ListElement *next;
};
struct PriorityQueue {
    ListElement *list;
};

bool pq_less_than(Node *n1, Node *n2) {
    if (n1->weight < n2->weight)
        return true;
    if (n1->weight > n2->weight)
        return false;
    return n1->symbol < n2->symbol;
}

PriorityQueue *pq_create(void) {

    PriorityQueue *pq = calloc(1, sizeof(PriorityQueue));

    return pq;
}

void pq_free(PriorityQueue **q) {

    free(*q);

    *q = NULL;
}

bool pq_is_empty(PriorityQueue *q) {

    return q->list == NULL;
}

bool pq_size_is_1(PriorityQueue *q) {

    return !(pq_is_empty(q)) && q->list->next == NULL;
}

void enqueue(PriorityQueue *q, Node *tree) {

    ListElement *e = calloc(1, sizeof(ListElement));

    e->tree = tree;

    if (pq_is_empty(q)) {

        q->list = e;

    } else if (pq_less_than(e->tree, q->list->tree)) {

        e->next = q->list;
        q->list = e;

    } else {
        ListElement *li = q->list;
        while (li->next != NULL && !pq_less_than(tree, li->next->tree)) {
            li = li->next;
        }

        e->next = li->next;
        li->next = e;
    }
}

bool dequeue(PriorityQueue *q, Node **tree) {
    if (pq_is_empty(q))
        return NULL;
    ListElement *head = q->list;

    q->list = head->next;

    *tree = head->tree;
    free(head);
    return true;
}

void pq_print(PriorityQueue *q) {
    assert(q != NULL);
    ListElement *e = q->list;
    int position = 1;
    while (e != NULL) {
        if (position++ == 1) {
            printf("=============================================\n");
        } else {
            printf("---------------------------------------------\n");
        }
        node_print_tree(e->tree, '<', 2);
        e = e->next;
    }
    printf("=============================================\n");
}
