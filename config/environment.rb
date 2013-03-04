# Load the rails application
require File.expand_path('../application', __FILE__)

# Initialize the rails application
Disicc::Application.initialize!

class Array
  def to_range
    array = self.compact.uniq.sort
    ranges = []
    if !array.empty?
      # Initialize the left and right endpoints of the range
      left, right = self.first, nil
      array.each do |obj|
        # If the right endpoint is set and obj is not equal to right's successor 
        # then we need to create a range.
        if right && obj != right.succ
          if right == left
            ranges << "#{right}"
          else
            ranges << "#{left}-#{right}"
          end
          left = obj
        end
        right = obj
      end
      if right == left
        ranges << "#{right}"
      else
        ranges << "#{left}-#{right}"
      end
    end
    ranges
  end
end
