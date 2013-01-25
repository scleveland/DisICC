class DisorderType 
  include DataMapper::Resource
  
  property :disorder_id, Serial
  property :disorder_type, String, :required => true
  property :deleted_at, ParanoidDateTime
  
end
