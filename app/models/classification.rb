class Classification
  include DataMapper::Resource
  
  property :class_id, Serial
  property :seq_id, Integer, :required => true
  property :class_name, String, :required => true
  
  # has n, :disorder_values
  # belongs_to :sequence
  
end
