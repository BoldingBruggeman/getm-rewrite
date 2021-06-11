class FortranObject:
    def get_arrays(self, names=(), halo=None, **kwargs):
        for name in names:
            data_ = self.get_array(name.encode('ascii'), **kwargs)

            # Store field with and without halo as class attribute
            setattr(self, name + '_', data_)
            setattr(self, name, data_[tuple([slice(halo, -halo)] * data_.ndim)])