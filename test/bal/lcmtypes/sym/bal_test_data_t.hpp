/** THIS IS AN AUTOMATICALLY GENERATED FILE.  DO NOT MODIFY
 * BY HAND!!
 *
 * Generated by lcm-gen
 *
 * From Source File: /Users/ryan/dev/symforce-zig-build/lcmtypes/test.lcm
 **/

#include <lcm/lcm_coretypes.h>

#ifndef __sym_bal_test_data_t_hpp__
#define __sym_bal_test_data_t_hpp__

#include <ostream>
#include <vector>
#include "lcmtypes/eigen_lcm/VectorXd.hpp"
#include "lcmtypes/sym/key_t.hpp"
#include "lcmtypes/sym/sparse_matrix_structure_t.hpp"

namespace sym
{

class bal_test_data_t
{
    public:
        std::vector< ::sym::key_t > optimized_keys;

        ::eigen_lcm::VectorXd residual;

        ::sym::sparse_matrix_structure_t hessian_structure;

        ::eigen_lcm::VectorXd hessian_data;

        ::eigen_lcm::VectorXd rhs;

        // solver stuff
        // hessian \ rhs
        ::eigen_lcm::VectorXd solution;

    public:
        bal_test_data_t() = default;

        /**
         * Member constructor
         */
        inline bal_test_data_t(
            const std::vector< ::sym::key_t >& optimized_keys_arg,
            const ::eigen_lcm::VectorXd& residual_arg,
            const ::sym::sparse_matrix_structure_t& hessian_structure_arg,
            const ::eigen_lcm::VectorXd& hessian_data_arg,
            const ::eigen_lcm::VectorXd& rhs_arg,
            const ::eigen_lcm::VectorXd& solution_arg
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

        using type_name_array_t = const char[16];

        inline static constexpr type_name_array_t* getTypeNameArrayPtr();

        /**
         * Returns "bal_test_data_t"
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
                if(fp->v == bal_test_data_t::getHash)
                    return 0;
            const __lcm_hash_ptr cp = { p, bal_test_data_t::getHash };

            uint64_t hash = 0x303af12e3b7d9479LL +
                ::sym::key_t::_computeHash(&cp) +
         ::eigen_lcm::VectorXd::_computeHash(&cp) +
         ::sym::sparse_matrix_structure_t::_computeHash(&cp) +
         ::eigen_lcm::VectorXd::_computeHash(&cp) +
         ::eigen_lcm::VectorXd::_computeHash(&cp) +
         ::eigen_lcm::VectorXd::_computeHash(&cp);

            return (hash<<1) + ((hash>>63)&1);
        }

        // Comparison operators.
        inline bool operator==(const bal_test_data_t& other) const;
        inline bool operator!=(const bal_test_data_t& other) const;

        // Ability to print to standard streams as well as the fmt library.
        friend std::ostream& operator<<(std::ostream& stream, const bal_test_data_t& obj) {
#if defined(SKYMARSHAL_PRINTING_ENABLED)
            stream << "bal_test_data_t(";
            stream << "optimized_keys=[";
            for (size_t i = 0; i < obj.optimized_keys.size(); ++i) {
                stream << obj.optimized_keys[i];
                if (i + 1 < obj.optimized_keys.size()) {
                    stream << ", ";
                }
            }
            stream << "]" << ", ";
            stream << "residual=" << obj.residual << ", ";
            stream << "hessian_structure=" << obj.hessian_structure << ", ";
            stream << "hessian_data=" << obj.hessian_data << ", ";
            stream << "rhs=" << obj.rhs << ", ";
            stream << "solution=" << obj.solution;
            stream << ")";
#else
            stream << "<FORMATTING DISABLED>";
#endif
            return stream;
        }
};

bal_test_data_t::bal_test_data_t(
    const std::vector< ::sym::key_t >& optimized_keys_arg,
    const ::eigen_lcm::VectorXd& residual_arg,
    const ::sym::sparse_matrix_structure_t& hessian_structure_arg,
    const ::eigen_lcm::VectorXd& hessian_data_arg,
    const ::eigen_lcm::VectorXd& rhs_arg,
    const ::eigen_lcm::VectorXd& solution_arg
) : optimized_keys(optimized_keys_arg),
    residual(residual_arg),
    hessian_structure(hessian_structure_arg),
    hessian_data(hessian_data_arg),
    rhs(rhs_arg),
    solution(solution_arg) {}

__lcm_buffer_size bal_test_data_t::encode(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;
    int64_t hash = (int64_t)getHash();

    tlen = __int64_t_encode_array(buf, offset + pos, maxlen - pos, &hash, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->_encodeNoHash(buf, offset + pos, maxlen - pos);
    if (tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size bal_test_data_t::decode(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
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

__lcm_buffer_size bal_test_data_t::getEncodedSize() const
{
    return 8 + _getEncodedSizeNoHash();
}

int64_t bal_test_data_t::getHash()
{
    static int64_t hash = _computeHash(NULL);
    return hash;
}

constexpr bal_test_data_t::type_name_array_t* bal_test_data_t::getTypeNameArrayPtr() {
    return &"bal_test_data_t";
}

constexpr const char* bal_test_data_t::getTypeName()
{
    return *bal_test_data_t::getTypeNameArrayPtr();
}

constexpr const char * bal_test_data_t::getPackageName()
{
    return "sym";
}

__lcm_buffer_size bal_test_data_t::_encodeNoHash(void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen) const
{
    __lcm_buffer_size pos = 0, tlen;

    int32_t v_num_optimized_keys = this->optimized_keys.size();
    tlen = __int32_t_encode_array(buf, offset + pos, maxlen - pos, &v_num_optimized_keys, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    for (__lcm_buffer_size a0 = 0; a0 < v_num_optimized_keys; a0++) {
        tlen = this->optimized_keys[a0]._encodeNoHash(buf, offset + pos, maxlen - pos);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    tlen = this->residual._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->hessian_structure._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->hessian_data._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->rhs._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->solution._encodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size bal_test_data_t::_decodeNoHash(const void *buf, __lcm_buffer_size offset, __lcm_buffer_size maxlen)
{
    __lcm_buffer_size pos = 0, tlen;

    int32_t v_num_optimized_keys;
    tlen = __int32_t_decode_array(buf, offset + pos, maxlen - pos, &v_num_optimized_keys, 1);
    if(tlen < 0) return tlen; else pos += tlen;

    this->optimized_keys.resize(v_num_optimized_keys);
    for (__lcm_buffer_size a0 = 0; a0 < v_num_optimized_keys; a0++) {
        tlen = this->optimized_keys[a0]._decodeNoHash(buf, offset + pos, maxlen - pos);
        if(tlen < 0) return tlen; else pos += tlen;
    }

    tlen = this->residual._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->hessian_structure._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->hessian_data._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->rhs._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    tlen = this->solution._decodeNoHash(buf, offset + pos, maxlen - pos);
    if(tlen < 0) return tlen; else pos += tlen;

    return pos;
}

__lcm_buffer_size bal_test_data_t::_getEncodedSizeNoHash() const
{
    __lcm_buffer_size enc_size = 0;
    enc_size += __int32_t_encoded_array_size(NULL, 1);
    for (__lcm_buffer_size a0 = 0; a0 < this->optimized_keys.size(); a0++) {
        enc_size += this->optimized_keys[a0]._getEncodedSizeNoHash();
    }
    enc_size += this->residual._getEncodedSizeNoHash();
    enc_size += this->hessian_structure._getEncodedSizeNoHash();
    enc_size += this->hessian_data._getEncodedSizeNoHash();
    enc_size += this->rhs._getEncodedSizeNoHash();
    enc_size += this->solution._getEncodedSizeNoHash();
    return enc_size;
}

bool bal_test_data_t::operator==(const bal_test_data_t& other) const {
  return (
          (optimized_keys==other.optimized_keys) && 
          (residual==other.residual) && 
          (hessian_structure==other.hessian_structure) && 
          (hessian_data==other.hessian_data) && 
          (rhs==other.rhs) && 
          (solution==other.solution));
}

bool bal_test_data_t::operator!=(const bal_test_data_t& other) const {
  return !(*this==other);
}

}

#endif