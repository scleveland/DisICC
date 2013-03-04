class Conservation
  include DataMapper::Resource
  
  property :id, Serial, :key=>true
  property :alignment_id, Integer, :required => true, :key=>true
  property :conservation_results, Text, :required => true
  property :deleted_at, ParanoidDateTime
  
  
end
