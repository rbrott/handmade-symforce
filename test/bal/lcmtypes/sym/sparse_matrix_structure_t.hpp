/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY
 * BY HAND!!
 *
 * Generated by lcm-gen
 *
 * From Source File: /Users/ryan/dev/symforce-zig-build/lcmtypes/symforce.lcm
 **/

#include <lcm/lcm_coretypes.h>

#ifndef __sym_sparse_matrix_structure_t_hpp__
#define __sym_sparse_matrix_structure_t_hpp__

#include <ostream>
#include <array>
#include "lcmtypes/eigen_lcm/VectorXi.hpp"

namespace sym
{

/**
 * The structure of a sparse matrix in CSC format, not including the numerical values
 * For a description of the format, see
 * https://en.wikipedia.org/wiki/Sparse_matrix#Compressed_sparse_column_(CSC_or_CCS)
 * In the comments below, assume an M x N matrix with nnz nonzeros
 */
class sparse_matrix_structure_t
{
    public:
        // The row for each nonzero entry in the matrix
        // size nnz
        ::eigen_lcm::VectorXi row_indices;

        // The index into row_indices (and the values) of the start of each column
        // size N
        ::eigen_lcm::VectorXi column_pointers;

        std::array< int64_t, 2 > shape;

    public:
        sparse_matrix_structure_t() = default;

        /**
         * Member constructor
         */
        inline sparse_matrix_structure_t(
            const ::eigen_lcm::VectorXi& row_indices_arg,
            const ::eigen_lcm::VectorXi& column_pointers_arg,
            const std::array< int64_t, 2 >& shape_arg
        );

        /**
         * Encode a message into binary form.
         *
         * @param buf The output buffer.
         * @param offset Encoding starts at thie byte offset into @p buf.
         * @param maxlen Maximum number of bytes to write.  This should generally be
         *  equal to getEncodedSize().
         * @return The number of bytes encoded, or <0 on error.
         */
        inline __lcm_buffer_size encode(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const;

        /**
         * Check how many bytes are required to encode this message.
         */
        inline __lcm_buffer_size getEncodedSize() const;

        /**
         * Decode a message from binary form into this instance.
         *
         * @param buf The buffer containing the encoded message.
         * @param offset The byte offset into @p buf where the encoded message starts.
         * @param maxlen The maximum number of bytes to read while decoding.
         * @return The number of bytes decoded, or <0 if an error occured.
         */
        inline __lcm_buffer_size decode(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen);

        /**
         * Retrieve the 64-bit fingerprint identifying the structure of the message.
         * Note that the fingerprint is the same for all instances of the same
         * message type, and is a fingerprint on the message type definition, not on
         * the message contents.
         */
        inline static int64_t getHash();

        using type_name_array_t = const char[26];

        inline static constexpr type_name_array_t* getTypeNameArrayPtr();

        /**
         * Returns "sparse_matrix_structure_t"
         */
        inline static constexpr const char* getTypeName();

        /**
         * Returns "sym"
         */
        inline static constexpr const char * getPackageName();

        // LCM support functions. Users should not call these
        inline __lcm_buffer_size _encodeNoHash(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const;
        inline __lcm_buffer_size _getEncodedSizeNoHash() const;
        inline __lcm_buffer_size _decodeNoHash(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen);
#if !defined(SKYDIO_DISABLE_LCM_NO_INLINE)
        __attribute__((noinline))
#endif
        static uint64_t _computeHash(const __lcm_hash_ptr *p)
        {
            const __lcm_hash_ptr *fp;
            for(fp = p; fp != NULL; fp = fp->parent)
                if(fp->v == sparse_matrix_structure_t::getHash)
                    return 0;
            const __lcm_hash_ptr cp = { p, sparse_matrix_structure_t::getHash };

            uint64_t hash = 0x5de55f243244b8acLL +
                ::eigen_lcm::VectorXi::_computeHash(&cp) +
         ::eigen_lcm::VectorXi::_computeHash(&cp);

            return (hash<<1) + ((hash>>63)&1);
        }

        // Comparison operators.
        inline bool operator==(const sparse_matrix_structure_t& other) const;
        inline bool operator!=(const sparse_matrix_structure_t& other) const;

        // Ability to print to standard streams as well as the fmt library.
        friend std::ostream& operator<<(std::ostream& stream, const sparse_matrix_structure_t& obj) {
#if defined(SKYMARSHAL_PRINTING_ENABLED)
            stream << "sparse_matrix_structure_t(";
            stream << "row_indices=" << obj.row_indices << ", ";
            stream << "column_pointers=" << obj.column_pointers << ", ";
            stream << "shape=[";
            for (size_t i = 0; i < obj.shape.size(); ++i) {
                stream << obj.shape[i];
                if (i + 1 < obj.shape.size()) {
                    stream << ", ";
                }
            }
            stream << "]";
            stream << ")";
#else
            stream << "<FORMATTING DISABLED>";
#endif
            return stream;
        }
};

sparse_matrix_structure_t::sparse_matrix_structure_t(
    const ::eigen_lcm::VectorXi& row_indices_arg,
    const ::eigen_lcm::VectorXi& column_pointers_arg,
    const std::array< int64_t, 2 >& shape_arg
) : row_indices(row_indices_arg),
    column_pointers(column_pointers_arg),
    shape(shape_arg) {}

__lcm_buffer_size sparse_matrix_structure_t::encode(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;
    int64_t hash = (int64_t)getHash();

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &hash, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->_encodeNoHash(buf, offset + pos, maxlen - pos);
    if (tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size sparse_matrix_structure_t::decode(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
{
    __lcm_buffer_size pos = 0, thislen;

    int64_t msg_hash;
    thislen = __int64_t_decode_array(buf, offset + pos, maxlen - pos, &msg_hash, 1);
    if (thislen < 0) return thislen; else pos += thislen;
    if (msg_hash != getHash()) return -1;

    thislen = this->_decodeNoHash(buf, offset + pos, maxlen - pos);
    if (thislen < 0) return thislen; else pos += thislen;

    return pos;
}

__lcm_buffer_size sparse_matrix_structure_t::getEncodedSize() const
{
    return 8 + _getEncodedSizeNoHash();
}

int64_t sparse_matrix_structure_t::getHash()
{
    static int64_t hash = _computeHash(NULL);
    return hash;
}

constexpr sparse_matrix_structure_t::type_name_array_t* sparse_matrix_structure_t::getTypeNameArrayPtr() {
    return &"sparse_matrix_structure_t";
}

constexpr const char* sparse_matrix_structure_t::getTypeName()
{
    return *sparse_matrix_structure_t::getTypeNameArrayPtr();
}

constexpr const char * sparse_matrix_structure_t::getPackageName()
{
    return "sym";
}

__lcm_buffer_size sparse_matrix_structure_t::_encodeNoHash(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;

    tlen = this->row_indices._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->column_pointers._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &this->shape[0], 2);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size sparse_matrix_structure_t::_decodeNoHash(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
{
    __lcm_buffer_size pos = 0, tlen;

    tlen = this->row_indices._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->column_pointers._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int64_t_decode_array(buf, offset + pos, maxlen - pos, &this->shape[0], 2);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size sparse_matrix_structure_t::_getEncodedSizeNoHash() const
{
    __lcm_buffer_size enc_size = 0;
    enc_size += this->row_indices._getEncodedSizeNoHash();
    enc_size += this->column_pointers._getEncodedSizeNoHash();
    enc_size += __int64_t_encoded_array_size(NULL, 2);
    return enc_size;
}

bool sparse_matrix_structure_t::operator==(const sparse_matrix_structure_t& other) const {
  return (
          (row_indices==other.row_indices) && 
          (column_pointers==other.column_pointers) && 
          (shape==other.shape));
}

bool sparse_matrix_structure_t::operator!=(const sparse_matrix_structure_t& other) const {
  return !(*this==other);
}

}

#endif