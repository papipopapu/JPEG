#include "image_compression.h"

bool is_little_endian() {
    uint16_t a = 1;
    return *(char *)&a == 1;
}

void OUTSTREAM_reset_bytes(OUTSTREAM *out) {
    if (out->written_bytes == (out->buffer_bytes-1)) { // there has been slip, and at we had alredy written up to the second to last byte
        fwrite(out->buffer, 1, out->buffer_bytes, out->file);
        memset(out->buffer, 0, out->buffer_bytes);
        out->written_bytes = -1;
    }   out->written_bits = 0;
        out->written_bytes++;
}
void INSTREAM_reset_bytes(INSTREAM *in) {
    if (in->read_bytes == (in->buffer_bytes-1)) { // there has been slip, and at we had alredy written up to the second to last byte
        in->buffer_bytes = fread(in->buffer, 1, in->buffer_bytes, in->file);
        if (in->buffer_bytes == 0) {in->eof = true;}
        in->read_bytes = -1;    
    }   in->read_bits = 0;
        in->read_bytes++;
        
}

int OUTSTREAM_push(OUTSTREAM *out, uint16_t data, int bits) {
    // bits are the number of bits to encode from least to most significant
    // printf("Pushing %d bits. Written bytes: %d, written bits: %d\n", bits, out->written_bytes, out->written_bits);
    if (bits < 1 || bits > 32) return 1;
    data &= (bits == 16) ? ~0 : (1 << bits) - 1; // remove trash data, but if whole int is needed, we cant bitshift by 0
    int slip = bits + out->written_bits - 8;  
    if (slip > 0) {
        out->buffer[out->written_bytes] |= data >> +slip;   
        OUTSTREAM_reset_bytes(out);
        OUTSTREAM_push(out, data, slip);
    } else if (slip < 0) {
        out->buffer[out->written_bytes] |= data << -slip;
        out->written_bits += bits; 
        // printf("< pushed %d data: ", bits); print_16bits(data << -slip); printf("\n");
    } else if (slip == 0) {
        out->buffer[out->written_bytes] |= data;
        // printf("== pushed %d data: ", bits); print_16bits(data); printf("\n");
        OUTSTREAM_reset_bytes(out);
    } 
    return 0;
}
int INSTREAM_pull_1bit(INSTREAM *in) {// returs 1 or 0 if successful, -1 if not

    if (in->eof) return -1;
    if (in->read_bits == 7) {
        int ret = in->buffer[in->read_bytes] & 1;
        INSTREAM_reset_bytes(in);
        return ret;
    } else {
        in->read_bits++;
        return (in->buffer[in->read_bytes] >> (8-in->read_bits)) & 1;//0000 0000
    }
}

int INSTREAM_pull(INSTREAM *in, uint16_t *data, int bits) { // pointer to uint, not array of uints

////////////////7 make this a dderivative from pull 1 bit dumb fuck retard
    // print buffer
    
    
    // please set data t 0 to avoid garbage data
    if (bits < 1 || bits > 16) return 1;
    if (in->eof) return -1;
    int i = 1;
    *data = INSTREAM_pull_1bit(in);
    while(i<bits) {
        *data <<= 1;
        *data |= INSTREAM_pull_1bit(in);
        i++;
    }
    return 0;
}
OUTSTREAM *new_OUTSTREAM(const char* filename, int buffer_bytes) {
    if (buffer_bytes < 1) return NULL;
    OUTSTREAM *out = malloc(sizeof(OUTSTREAM));
    out->filename = filename;
    out->file = fopen(filename, "wb");
    out->written_bits = 0;
    out->buffer_bytes = buffer_bytes;
    out->buffer = calloc(buffer_bytes, 1); 
    //if (fread(out->buffer, 1, buffer_bytes, out->file)==0) {fclose(out->file); free(out); return NULL;}
    return out;
}
INSTREAM *new_INSTREAM(const char* filename, int buffer_bytes) {
    if (buffer_bytes < 1) return NULL;
    INSTREAM *in = malloc(sizeof(INSTREAM));
    in->filename = filename;
    in->file = fopen(filename, "rb");
    in->read_bits = 0;
    in->read_bytes = 0;
    in->eof = false;
    in->buffer_bytes = buffer_bytes;
    in->buffer = calloc(buffer_bytes, 1); 
    in->buffer_bytes = fread(in->buffer, 1, buffer_bytes, in->file);
    if (in->buffer_bytes == 0) {fclose(in->file); free(in->buffer); free(in); return NULL;}
    return in;
}

int delete_OUTSTREAM(OUTSTREAM *out) {
    if (out == NULL) return -1;
    if (out->written_bits != 0 || out->written_bytes != 0) fwrite(out->buffer, 1, out->written_bytes+1, out->file);
    // printf("End writing %d bytes and %d bits. Buffer:", out->written_bytes, out->written_bits); 
    // print at least 1 uint from buffer

// flush buffer if something has been written
    fclose(out->file);
    free(out->buffer);
    free(out);
    return 0;
}
int delete_INSTREAM(INSTREAM *in) {
    if (!in) return -1;
    fclose(in->file);
    free(in->buffer);
    free(in);
    return 0;
}