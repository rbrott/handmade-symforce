/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY
 * BY HAND!!
 *
 * Generated by lcm-gen
 *
 * From Source File: /Users/ryan/dev/symforce-zig-build/lcmtypes/symforce.lcm
 **/

#include <lcm/lcm_coretypes.h>

#ifndef __sym_linearization_sparse_factor_helper_t_hpp__
#define __sym_linearization_sparse_factor_helper_t_hpp__

#include <ostream>
#include <vector>
#include "lcmtypes/sym/linearization_offsets_t.hpp"

namespace sym
{

class linearization_sparse_factor_helper_t
{
    public:
        // Total residual dimension of the factor
        int32_t residual_dim;

        // Offset of this factor's residual slice within the whole problem residual
        int32_t combined_residual_offset;

        std::vector< ::sym::linearization_offsets_t > key_helpers;

        std::vector< int32_t > jacobian_index_map;

        std::vector< int32_t > hessian_index_map;

    public:
        linearization_sparse_factor_helper_t() = default;

        /**
         * Member constructor
         */
        inline linearization_sparse_factor_helper_t(
            const int32_t& residual_dim_arg,
            const int32_t& combined_residual_offset_arg,
            const std::vector< ::sym::linearization_offsets_t >& key_helpers_arg,
            const std::vector< int32_t >& jacobian_index_map_arg,
            const std::vector< int32_t >& hessian_index_map_arg
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

        using type_name_array_t = const char[37];

        inline static constexpr type_name_array_t* getTypeNameArrayPtr();

        /**
         * Returns "linearization_sparse_factor_helper_t"
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
                if(fp->v == linearization_sparse_factor_helper_t::getHash)
                    return 0;
            const __lcm_hash_ptr cp = { p, linearization_sparse_factor_helper_t::getHash };

            uint64_t hash = 0xec50286eabbb5f0bLL +
                ::sym::linearization_offsets_t::_computeHash(&cp);

            return (hash<<1) + ((hash>>63)&1);
        }

        // Comparison operators.
        inline bool operator==(const linearization_sparse_factor_helper_t& other) const;
        inline bool operator!=(const linearization_sparse_factor_helper_t& other) const;

        // Ability to print to standard streams as well as the fmt library.
        friend std::ostream& operator<<(std::ostream& stream, const linearization_sparse_factor_helper_t& obj) {
#if defined(SKYMARSHAL_PRINTING_ENABLED)
            stream << "linearization_sparse_factor_helper_t(";
            stream << "residual_dim=" << obj.residual_dim << ", ";
            stream << "combined_residual_offset=" << obj.combined_residual_offset << ", ";
            stream << "key_helpers=[";
            for (size_t i = 0; i < obj.key_helpers.size(); ++i) {
                stream << obj.key_helpers[i];
                if (i + 1 < obj.key_helpers.size()) {
                    stream << ", ";
                }
            }
            stream << "]" << ", ";
            stream << "jacobian_index_map=[";
            for (size_t i = 0; i < obj.jacobian_index_map.size(); ++i) {
                stream << obj.jacobian_index_map[i];
                if (i + 1 < obj.jacobian_index_map.size()) {
                    stream << ", ";
                }
            }
            stream << "]" << ", ";
            stream << "hessian_index_map=[";
            for (size_t i = 0; i < obj.hessian_index_map.size(); ++i) {
                stream << obj.hessian_index_map[i];
                if (i + 1 < obj.hessian_index_map.size()) {
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

linearization_sparse_factor_helper_t::linearization_sparse_factor_helper_t(
    const int32_t& residual_dim_arg,
    const int32_t& combined_residual_offset_arg,
    const std::vector< ::sym::linearization_offsets_t >& key_helpers_arg,
    const std::vector< int32_t >& jacobian_index_map_arg,
    const std::vector< int32_t >& hessian_index_map_arg
) : residual_dim(residual_dim_arg),
    combined_residual_offset(combined_residual_offset_arg),
    key_helpers(key_helpers_arg),
    jacobian_index_map(jacobian_index_map_arg),
    hessian_index_map(hessian_index_map_arg) {}

__lcm_buffer_size linearization_sparse_factor_helper_t::encode(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;
    int64_t hash = (int64_t)getHash();

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &hash, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->_encodeNoHash(buf, offset + pos, maxlen - pos);
    if (tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size linearization_sparse_factor_helper_t::decode(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
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

__lcm_buffer_size linearization_sparse_factor_helper_t::getEncodedSize() const
{
    return 8 + _getEncodedSizeNoHash();
}

int64_t linearization_sparse_factor_helper_t::getHash()
{
    static int64_t hash = _computeHash(NULL);
    return hash;
}

constexpr linearization_sparse_factor_helper_t::type_name_array_t* linearization_sparse_factor_helper_t::getTypeNameArrayPtr() {
    return &"linearization_sparse_factor_helper_t";
}

constexpr const char* linearization_sparse_factor_helper_t::getTypeName()
{
    return *linearization_sparse_factor_helper_t::getTypeNameArrayPtr();
}

constexpr const char * linearization_sparse_factor_helper_t::getPackageName()
{
    return "sym";
}

__lcm_buffer_size linearization_sparse_factor_helper_t::_encodeNoHash(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->residual_dim, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->combined_residual_offset, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    int32_t v_num_key_helpers = this->key_helpers.size();
    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &v_num_key_helpers, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    for (__lcm_buffer_size a0 = 0; a0 < v_num_key_helpers; a0++) {
        tlen = this->key_helpers[a0]._encodeNoHash(buf, offset + pos, maxlen - pos);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    int32_t v_num_jacobian_index_map = this->jacobian_index_map.size();
    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &v_num_jacobian_index_map, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    if(v_num_jacobian_index_map > 0) {
        tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->jacobian_index_map[0], v_num_jacobian_index_map);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    int32_t v_num_hessian_index_map = this->hessian_index_map.size();
    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &v_num_hessian_index_map, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    if(v_num_hessian_index_map > 0) {
        tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &this->hessian_index_map[0], v_num_hessian_index_map);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    return pos;
}

__lcm_buffer_size linearization_sparse_factor_helper_t::_decodeNoHash(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
{
    __lcm_buffer_size pos = 0, tlen;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->residual_dim, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->combined_residual_offset, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    int32_t v_num_key_helpers;
    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &v_num_key_helpers, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    this->key_helpers.resize(v_num_key_helpers);
    for (__lcm_buffer_size a0 = 0; a0 < v_num_key_helpers; a0++) {
        tlen = this->key_helpers[a0]._decodeNoHash(buf, offset + pos, maxlen - pos);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    int32_t v_num_jacobian_index_map;
    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &v_num_jacobian_index_map, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    this->jacobian_index_map.resize(v_num_jacobian_index_map);
    if(v_num_jacobian_index_map) {
        tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->jacobian_index_map[0], v_num_jacobian_index_map);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    int32_t v_num_hessian_index_map;
    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &v_num_hessian_index_map, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    this->hessian_index_map.resize(v_num_hessian_index_map);
    if(v_num_hessian_index_map) {
        tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &this->hessian_index_map[0], v_num_hessian_index_map);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    return pos;
}

__lcm_buffer_size linearization_sparse_factor_helper_t::_getEncodedSizeNoHash() const
{
    __lcm_buffer_size enc_size = 0;
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    for (__lcm_buffer_size a0 = 0; a0 < this->key_helpers.size(); a0++) {
        enc_size += this->key_helpers[a0]._getEncodedSizeNoHash();
    }
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    enc_size += __int32_t_encoded_array_size(NULL, this->jacobian_index_map.size());
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    enc_size += __int32_t_encoded_array_size(NULL, this->hessian_index_map.size());
    return enc_size;
}

bool linearization_sparse_factor_helper_t::operator==(const linearization_sparse_factor_helper_t& other) const {
  return (
          (residual_dim==other.residual_dim) && 
          (combined_residual_offset==other.combined_residual_offset) && 
          (key_helpers==other.key_helpers) && 
          (jacobian_index_map==other.jacobian_index_map) && 
          (hessian_index_map==other.hessian_index_map));
}

bool linearization_sparse_factor_helper_t::operator!=(const linearization_sparse_factor_helper_t& other) const {
  return !(*this==other);
}

}

#endif