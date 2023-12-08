

#include "bitwriter.h"
#include "pq.h"

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#define OPTIONS "i:o:h"
typedef struct Code {
    uint64_t code;
    uint8_t code_length;
} Code;

uint64_t fill_histogram(Buffer *inbuf, double *histogram) {
    uint8_t byte;

    for (int i = 0; i < 256; i++) {
        histogram[i] = 0.0;
    }

    // Initialize the file size counter
    uint64_t filesize = 0;
    histogram[0x00]++;
    histogram[0xff]++;
    // Read bytes from the input buffer and update the histogram
    while (read_uint8(inbuf, &byte)) {
        // read_uint8(inbuf, &byte);
        histogram[byte]++;
        filesize++;
    }

    // Apply the important hack to ensure at least two non-zero values in the histogram
    // histogram[0x00]++;
    // histogram[0xff]++;

    return filesize;
}

Node *create_tree(double *histogram, uint16_t *num_leaves) {

    PriorityQueue *pq = pq_create();
    // init_priority_queue(&pq);
    uint16_t x;
    for (int i = 0; i < 256; i++) {
        if (histogram[i] != 0.0) {
            Node *node = node_create(i, histogram[i]);
            enqueue(pq, node);
            (*num_leaves)++;
        }
    }

    // Step 2: Run the Huffman Coding algorithm
    while (!pq_size_is_1(pq) && !pq_is_empty(pq)) {
        Node *left;
        Node *right;
        dequeue(pq, &left);
        dequeue(pq, &right);

        Node *new_node = node_create(0x00, left->weight + right->weight);
        new_node->left = left;
        new_node->right = right;

        enqueue(pq, new_node);
    }

    // Step 3: Dequeue the only entry remaining and return it as the Huffman Tree
    //    num_leaves = &x;
    Node *huffman_tree;

    dequeue(pq, &huffman_tree);
    pq_free(&pq);

    return huffman_tree;
}

// probably correct, check for the casting thing in front of the one
void fill_code_table(Code *code_table, Node *node, uint64_t code, uint8_t code_length) {
    //	code_table = calloc(256, sizeof(Code));
    if (node->left != NULL && node->right != NULL) {
        // Internal node
        fill_code_table(code_table, node->left, code, code_length + 1);

        code |= (1 << code_length);

        fill_code_table(code_table, node->right, code, code_length + 1);
    } else {
        // Leaf node: store the Huffman code
        code_table[node->symbol].code = code;
        code_table[node->symbol].code_length = code_length;
    }

    //	free(code_table);
}

void huff_write_tree(BitWriter *outbuf, Node *node) {
    if (node->left != NULL && node->right != NULL) {
        // Internal node
        huff_write_tree(outbuf, node->left);
        huff_write_tree(outbuf, node->right);
        bit_write_bit(outbuf, 0);
    } else {
        // Leaf node
        bit_write_bit(outbuf, 1);
        bit_write_uint8(outbuf, node->symbol);
    }
}

void node_clear(Node *n) {

    if (n->left) {
        node_clear(n->left);
    }

    if (n->right) {

        node_clear(n->right);
    }

    free(n);
}

void huff_compress_file(BitWriter *outbuf, Buffer *inbuf, uint32_t filesize, uint16_t num_leaves,
    Node *code_tree, Code *code_table) {

    uint8_t byte;

    bit_write_uint8(outbuf, 'H');
    bit_write_uint8(outbuf, 'C');
    bit_write_uint32(outbuf, filesize);
    bit_write_uint16(outbuf, num_leaves);

    // Write the code tree
    huff_write_tree(outbuf, code_tree);

    // Read the input file and write the Huffman codes
    while (read_uint8(inbuf, &byte)) {
        //        read_uint8(inbuf, &byte);
        //
        uint64_t code = code_table[byte].code;
        uint8_t code_length = code_table[byte].code_length;

        for (int i = 0; i <= code_length - 1; i++) {

            bit_write_bit(outbuf, code & 1);
            code >>= 1;
        }
    }
}

int main(int argc, char *argv[]) {

    //  FILE *infile = stdin;
    // FILE *outfile = stdout;

    char *file_in;
    char *file_out;

    //int i = 0;
    //int o = 0;
    int h = 0;

    int opt = 0;

    while ((opt = getopt(argc, argv, OPTIONS)) != -1) {

        switch (opt) {

        case 'i':
            file_in = optarg;
            if (file_in == NULL) {

                printf("error!");
                exit(1);
            }
            //    i = 1;
            break;

        case 'o':
            file_out = optarg;
            //   if (file_out == NULL) {
            //     printf("error!");
            //   exit(1);
            // }
            //printf("error!");
            // exit(1);
            //  o = 1;
            break;

        case 'h':
            h = 1;
            printf("i: set the name of the file to be read from.\no: set the name of the file to "
                   "be written to.\n h: print a help message.");
        }
    }

    if (h != 1) {

        Buffer *b = read_open(file_in);
        // Code *c = calloc(256, sizeof(Code));

        double histogram[256];
        // BitWriter *bw = bit_write_open(b);
        uint32_t filesize = (uint32_t) fill_histogram(b, histogram);
        uint16_t num_leaves = 0;

        Node *n = create_tree(histogram, &num_leaves);
        Code *c = calloc(256, sizeof(Code));

        // create_tree(histogram, &num_leaves);

        fill_code_table(c, n, 0, 0);
        read_close(&b);

        //open outfile into buffer
        b = read_open(file_in);
        BitWriter *out = bit_write_open(file_out);
        //huff compress with outfile and infile buffer
        huff_compress_file(out, b, filesize, num_leaves, n, c);
        //and then write_close and read_close()

        node_clear(n);
        free(c);

        read_close(&b);
        // write_close(&b1);
        // node_clear(n);
        bit_write_close(&out);
        // free(c);
        //        free(n);
    }
    return 0;
}
