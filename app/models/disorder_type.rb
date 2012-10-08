class DisorderType 
  include DataMapper::Resource
  
  property :disorder_id, Serial
  property :disorder_type, String, :required => true
  
end
